classdef ConConSmooth < matlab.mixin.SetGet
    %ConConSmooth: A class to implement a collection of fitted ConCon object
    
    properties
        Basis
        CoefMat
        reg_param
    end
    
    methods
        function obj = ConConSmooth(Basis,varargin)
            %ConConSmooth Construct an instance of this class
            % Input:
            %     Basis: ConConBasis object 
            %     varargin optional arguments
            %     CoefMat: N x K matrix of coefficient representations for
            %               Basis 
            %     reg_param: (float), speficying regulation strength
            
            p = inputParser;
            addParameter(p, 'CoefMat', [], @isnumeric)
            addParameter(p, 'reg_param', [], @isnumeric)
            parse(p, varargin{:});
            params = p.Results;
            
            obj.Basis = Basis;
            obj.CoefMat = params.CoefMat;
            obj.reg_param = params.reg_param;
        end
        
        function [CMat] = smooth(obj, Y, X0, X1)
            %smooth Compute continous representation from discrete
            %connectivity
            % Input:
            %      Y: length N cell array of sparse matrices (nverts,
            %      nverts)
            %      X0, X1: (nverts, 3) (euclidean) coordinates of the high-res
            
            N = length(Y);
            nv0 = size(X0,1);
            nv1 = size(X1,1);
            X = [X0;X1];
            ix = [zeros(1,nv0), ones(1,nv1)];
            Psi = obj.Basis.Evaluate(X, ix);
            K = get(obj.Basis, 'Rank');
            
            Psi_DM = zeros((nv0+nv1)*(nv0+nv1+1)/2, K);
            tri_ix = tril(true([nv0+nv1,nv0+nv1]));
            for k=1:K
                Psi_k = Psi(:,k)*Psi(:,k)';
                Psi_DM(:,k) = Psi_k(tri_ix);
            end 
            
            CMat = zeros(N,K);
            for i=1:N
                Y_i = Y{i};
                Y_i_vec = Y_i(tri_ix);
                CMat(i,:) = lsqr(Psi_DM, Y_i_vec, 1e-10, 20);
            end
            %obj.CoefMat = CMat;
        end     
    end
end

