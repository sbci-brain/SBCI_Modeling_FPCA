classdef ConConBasis < matlab.mixin.SetGet
    %ConConBasis: A class to implement a ConCon Basis Object 
    
    properties
        Basis0 %Marginal basis for a copy of S2 
        Basis1 %Marginal basis for a copy of S2 
        C0 % Cofficient matrix 0
        C1 % Coefficient matrix 1 
        Rank % Global Rank
    end
    
    methods
        function obj = ConConBasis(Basis0,Basis1,varargin)
            %ConConBasis Construct an instance of this class
            % Input:
            %     Basis0, Basis1: Basis object 
            %     varargin optional arguments
            %     K: rank of basis 
            
            p = inputParser;
            addParameter(p, 'K', 1, @isnumeric)
            parse(p, varargin{:});
            params = p.Results;
            
            obj.Rank = params.K;
            
            obj.Basis0 = Basis0;
            obj.Basis1 = Basis1;
            
        end
        function obj = set.C0(obj,value)
            obj.C0 = value;
        end 
        function obj = set.C1(obj,value)
            obj.C1 = value;
        end 
        function Psi = Evaluate(obj, X, ix)
            %Evaluate the ConCon object over a set of points
            %
            % Input:
            %   X: (Npoints, 3) (euclidean) coordinates of endpoints
            %   ix: (Npoints, 1) index into which copy of S2 (0,1), which
            %   must be ordered 
            %
            % Output:
            %    Psi: (Npoints, K) basis evaluation matrix 
            
            X0 = X(ix==0,:); X1 = X(ix==1,:);
            Phi0 = obj.Basis0.Evaluate(X0);
            Phi1 = obj.Basis1.Evaluate(X1);
            Psi0 = Phi0*obj.C0;
            Psi1 = Phi1*obj.C1;
            Psi = [Psi0;Psi1];
        end
        
        function GradPsi = MargGradient(obj, X, ix)
            %Construct the marginal gradient tensor of the ConCon basis 
            %
            % Input:
            %   X: (Npoints, 3) (euclidean) coordinates of endpoints
            %   ix: (Npoints, 1) index into which copy of S2 (0,1), which
            %   must be ordered 
            %
            % Output:
            %   GradPsi: (K, Npoints, 3) marginal gradient tensor field 
            X0 = X(ix==0,:); X1 = X(ix==1,:);
            GradPhi0 = tensor(obj.Basis0.EvalGrad(X0));
            GradPhi1 = tensor(obj.Basis1.EvalGrad(X1));
            GradPsi0 = double(ttm(GradPhi0, obj.C0', 1));
            GradPsi1 = double(ttm(GradPhi1, obj.C1', 1));
            GradPsi = [GradPsi0; GradPsi1];
        end 
        
        function GradConCon = JointGradient(obj, X1, ix1, X2, ix2)
            %Construct the gradient tensor of the ConCon basis over the
            %symmetric space \Omega x \Omega at each point in point-cloud
            % Input:
            %   X1: (Npoints, 3) (euclidean) coordinates of endpoint 1
            %   ix: (Npoints, 1) index into which copy of S2 (0,1), which
            %      must be ordered 
            %   X2: (Npoints, 3) (euclidean) coordinates of endpoint 2
            %   ix2: (Npoints, 1) index into which copy of S2 (0,1), which
            %      must be ordered 
            
            GradConCon = zeros(K, size(X1,1), 6);

            MGradPsi1 = obj.MargGradient(X1, ix1);
            Psi1 = obj.Evaluate(X1, ix1);
            MGradPsi2 = obj.MargGradient(X2, ix2);
            Psi2 = obj.Evaluate(X2, ix2);
             
            GradConCon(:,:,1) = MGradPsi1(:,:,1).*Psi2';
            GradConCon(:,:,2) = MGradPsi1(:,:,2).*Psi2';
            GradConCon(:,:,3) = MGradPsi1(:,:,3).*Psi2';
            
            GradConCon(:,:,4) = MGradPsi2(:,:,1).*Psi1';
            GradConCon(:,:,5) = MGradPsi2(:,:,2).*Psi1';
            GradConCon(:,:,6) = MGradPsi2(:,:,3).*Psi1';
        end 
        
        function [obj, Smat, scales, objective_function, times, residual_norms] = Fit(obj,Y,X0,X1,varargin)
            %Fit the coefficients to the observed data
            %
            % Input:
            %      Y: length N cell array of sparse matrices (nverts,
            %      nverts)
            %      X0, X1: (nverts, 3) (euclidean) coordinates of the high-res
            %      grid
            %      varargin optional arguments
            %      alpha_1: (float) smoothness penalty
            %      alpha_2: (int) local support constraint 
            %      MAXITER_OUTER: (int) max # of outer loop iterations 
            %      MAXITER_INNER: (int) max # of inner loop iterations 
            %      TOL_OUTER: (float) tolerance for iteration termination
            %      of outer loop
            %      TOL_INNER: (float) tolerance for iteration termination
            %      of inner loop
            %
            % Output:
            %
            
            % parse optional arguments    
            p = inputParser;
            addParameter(p, 'alpha_1', 1e-10, @isnumeric)
            addParameter(p, 'alpha_2', Inf, @isnumeric)
            addParameter(p, 'MAXITER_OUTER', 30, @isnumeric)
            addParameter(p, 'MAXITER_INNER', 30, @isnumeric)
            addParameter(p, 'TOL_OUTER', 1e-3, @isnumeric)
            addParameter(p, 'TOL_INNER', 1e-3, @isnumeric)
            parse(p, varargin{:});
            params = p.Results;
            
            % set global parameters 
            alpha_1 = params.alpha_1;
            alpha_2 = params.alpha_2;
            MAXITER_OUTER = params.MAXITER_OUTER;
            MAXITER_INNER = params.MAXITER_INNER;
            TOL_OUTER = params.TOL_OUTER;
            TOL_INNER = params.TOL_INNER;
            m0 = get(obj.Basis0, 'M');
            m1 = get(obj.Basis1, 'M');
            K = obj.Rank;
            nv0 = length(X0);
            nv1 = length(X1);
            N = length(Y);

            % constuct basis evaluation matrix and perform data transformation 
            Phi0 = obj.Basis0.Evaluate(X0);
            Phi1 = obj.Basis1.Evaluate(X1);

            [U0, S0, V0] = svd(Phi0, 'econ');
            [U1, S1, V1] = svd(Phi1, 'econ');
            U = [U0, zeros(size(U0, 1), size(U1,2)); 
                zeros(size(U1, 1), size(U0,2)), U1];
            
            U_sparse = sparse(U);
            
            % Perform data transformation 
            G = zeros(m0+m1, m0+m1, N);
            for i=1:N
                Y_i = Y{i};
                G(:,:,i) = full(U_sparse'*Y_i*U_sparse);
            end
            Gtensor = tensor(G);
            
            % build inner product matrices 
            obj.Basis0.InnerProduct()
            obj.Basis1.InnerProduct()
            obj.Basis0.QVInnerProduct()
            obj.Basis1.QVInnerProduct()
            
            % stack for inner prod on \Omega X \Omega
            J0 = get(obj.Basis0, 'J');
            J1 = get(obj.Basis1, 'J');
            R0 = get(obj.Basis0, 'R');
            R1 = get(obj.Basis1, 'R');
            J_xi =  [J0, zeros(m0, m1);
                    zeros(m1, m0), J1];
            R_xi = [R0, zeros(m0, m1);
                    zeros(m1, m0), R1];

            % transform inner product matrix             
            J_xi_tilde = [diag(1./diag(S0))*V0'*J0*V0*diag(1./diag(S0)), zeros(m0, m1);
                             zeros(m1, m0), diag(1./diag(S1))*V1'*J1*V1*diag(1./diag(S1))];
            R_xi_tilde = [diag(1./diag(S0))*V0'*R0*V0*diag(1./diag(S0)), zeros(m0, m1);
                         zeros(m1, m0), diag(1./diag(S1))*V1'*R1*V1*diag(1./diag(S1))];
                     
            % compute transformation matrices 
            A = blkdiag(V0*diag(1./diag(S0)),V1*diag(1./diag(S1)));
            Ainv = blkdiag(S0*V0',S1*V1');
            
            % global parameters for algorithm 
            m = m0+m1;

            % Pre-allocate storage 
            Ctildes = zeros(m,K);
            C = zeros(m,K);
            Smat = zeros(N,K);
            scales = zeros(K,1);
            times = zeros(K, MAXITER_OUTER, 2);
            objective_function = zeros(K, MAXITER_OUTER);
            residual_norms = zeros(K,1);
            
            % initialize tensor residual
            Ghat = Gtensor; 
            for k=1:K
            % Create projection matrix 
            if k == 1
                P = eye(m);
            else
                Ctilde_k = Ctildes(:,1:(k-1));
                Norm_mat = inv(Ctilde_k'*J_xi_tilde*Ctilde_k);
                P = eye(m) - Ctilde_k*Norm_mat*Ctilde_k'*J_xi_tilde;
            end 
            % Initialize factors 
            Gc = tenmat(Ghat,1);
            Gc2 = double(Gc * Gc');
            v0_c = normrnd(0,1, [size(Gc2,1), 1]);
            v0_c = v0_c/norm(v0_c);
            Vtraj_c0 = power_iterations(Gc2, v0_c, MAXITER_INNER, TOL_INNER);
            ctilde_k = Vtraj_c0(:,end);
            ctilde_k = ctilde_k/norm(ctilde_k);
            ghat_c0 = double(ttv(Ghat, ctilde_k, 1))'*ctilde_k;
            s_k = ghat_c0/norm(ghat_c0);
            % Evaluate objective 
            objective_function(k,1) = ctilde_k' * P *(double(ttv(Ghat, s_k, 3)) - alpha_1*R_xi_tilde) * P' * ctilde_k;
            i = 1; delta_residual_norm = inf;
            while (i < MAXITER_OUTER) && (delta_residual_norm > TOL_OUTER)
                % Update s 
                tic
                ghat_c = double(ttv(Ghat, ctilde_k, 1))'*ctilde_k;
                times(k, i, 1) = toc;
                s_k = ghat_c/norm(ghat_c);
                % Update c 
                Ghat_s = double(ttv(Ghat, s_k, 3));
                if isinf(alpha_2)
                    M_k = P *(Ghat_s - alpha_1*R_xi_tilde) * P';
                    M_k = (M_k + M_k')/2;
                    tic
                    [ctilde_k, d_] = eigs(M_k, 1);
                    times(k, i, 2) = toc;
                    c_k = A*ctilde_k/norm(A*ctilde_k);
                else
                    M_k = Ainv' * P *(Ghat_s - alpha_1*R_xi_tilde) * P' * Ainv;
                    v0 = A*ctilde_k;
                    tic
                    Vtrajectory = generalized_power(M_k, J_xi, v0, alpha_2, MAXITER_INNER, TOL_INNER);
                    times(k, i, 2) = toc;
                    ctilde_k = Ainv*Vtrajectory(:,end)/norm(Ainv*Vtrajectory(:,end));
                    c_k = Vtrajectory(:,end);
                end 
                % Compute objective function 
                objective_function(k,i+1) = ctilde_k' * P *(Ghat_s - alpha_1*R_xi_tilde) * P' * ctilde_k;
                delta_residual_norm = abs((objective_function(k,i+1) - objective_function(k,i))/objective_function(k,1));
                i = i + 1;
            end 
            % Compute scale and deflate the residual tensor 
            gamma = ctilde_k'*double(ttv(Ghat, s_k, 3))*ctilde_k;
            Ghat = Ghat - full(ktensor(gamma, ctilde_k, ctilde_k, s_k));
            % Record factors 
            C(:,k) = c_k;
            Ctildes(:,k) = ctilde_k;
            Smat(:,k) = s_k;
            scales(k) = gamma;
            % Compute the proportion of variance explained 
            residual_norms(k) = norm(Gtensor - Ghat)/norm(Gtensor);
            end 
           obj.C0 = C(1:m0,:);
           obj.C1 = C(m0+1:end,:);
        end 
    end
end

