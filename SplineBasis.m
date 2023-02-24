classdef SplineBasis < matlab.mixin.SetGet
    %SplineBasis: A class to implement a spherical spline basis 
    %   Detailed explanation goes here
    
    properties
        Verts
        TRI
        E 
        IE
        M
        D
        J
        R
    end
    
    methods
        function obj = SplineBasis(Verts)
            %Basis Construct an instance of this class
            % Input:
            %     Verts: (M, 3) (euclidean) coordinates of vertices, i.e.
            %     nodes of basis functions 
            [v1,v2,v3] = sdelaunay(Verts(:,1),Verts(:,2),Verts(:,3));
            TRI = [v1,v2,v3];
            [v1,v2,v3,e1,e2,e3,ie1,ie2,tril,trir] =  slists(Verts(:,1),Verts(:,2),Verts(:,3),TRI);
            
            obj.Verts = Verts;
            obj.TRI = TRI;
            obj.E = [e1,e2,e3];
            obj.IE = [ie1,ie2];
            obj.M = length(Verts(:,1));
            obj.D = 1; 
        end
        
        function Phi = Evaluate(obj,X)
            %Evaluate basis object evaluations at points in X 
            % Input:
            %       X: (Npoints, 3) (euclidean) coordinates on the sphere
            %      for evaluation 
            % Output:
            %       Phi: (Npoints, M) basis evaluation matrix                   
            x = obj.Verts(:,1);
            y = obj.Verts(:,2);
            z = obj.Verts(:,3);
            
            v1 = obj.TRI(:,1);
            v2 = obj.TRI(:,2);
            v3 = obj.TRI(:,3);
            
            e1 = obj.E(:,1);
            e2 = obj.E(:,2);
            e3 = obj.E(:,3);
            
            ie1 = obj.IE(:,1);
            
            m = obj.M;
            d = obj.D;
            
            nv = length(v1); Npoints = size(X,1); 
            Phi = zeros(Npoints, m);
            for np = 1:Npoints
                xp = X(np, 1); yp = X(np, 2); zp = X(np, 3);
                % get spherical triangle containing point 
                [it,b1,b2,b3] = sfindtri(x,y,z,v1,v2,v3,xp,yp,zp);
                % get basis functions assocaited with triangle 
                inds = getindex(d,it,m,v1,v2,v3,e1,e2,e3,ie1); 
                % evaluate each basis function at point 
                cb1 = double(1:m == inds(1));
                co1 = getco(d,it,nv,v1,v2,v3,e1,e2,e3,ie1,cb1); 
                val1 = decast(d,b1,b2,b3,co1);

                cb2 = double(1:m == inds(2));
                co2 = getco(d,it,nv,v1,v2,v3,e1,e2,e3,ie1,cb2); 
                val2 = decast(d,b1,b2,b3,co2);

                cb3 = double(1:m == inds(3));
                co3 = getco(d,it,nv,v1,v2,v3,e1,e2,e3,ie1,cb3); 
                val3 = decast(d,b1,b2,b3,co3);

                %store evaluations 
                Phi(np, inds(1))=val1; Phi(np,inds(2))=val2; Phi(np,inds(3))=val3;
            end 
        end
        
        function GradPhi = EvalGrad(obj, X)
            %Evaluate the gradient of the basis object at points in X 
            % Input:
            %       X: (Npoints, 3) (euclidean) coordinates on the sphere
            % 
            % Output:
            %       GradPhi: (M, Npoints, 3) basis gradient tensor
            
            x = obj.Verts(:,1);
            y = obj.Verts(:,2);
            z = obj.Verts(:,3);
            
            v1 = obj.TRI(:,1);
            v2 = obj.TRI(:,2);
            v3 = obj.TRI(:,3);
            
            e1 = obj.E(:,1);
            e2 = obj.E(:,2);
            e3 = obj.E(:,3);
            
            ie1 = obj.IE(:,1);
            
            m = obj.M;
            d = obj.D;
            
            
            nv = length(v1); Npoints = size(X,1);
            GradPhi = zeros(m,Npoints,3);
            for nq=1:Npoints
                % get triangle quadrature point is within 
                [it,b1,b2,b3] = sfindtri(x,y,z,v1,v2,v3,X(nq,1),X(nq,2),X(nq,3));
                % get basis functions assocaited with triangle 
                cb1 = double(1:m == v1(it)); cb2 = double(1:m == v2(it)); cb3 = double(1:m == v3(it));
                % get local coordiantes
                co1 = getco(d,it,nv,v1,v2,v3,e1,e2,e3,ie1,cb1); 
                co2 = getco(d,it,nv,v1,v2,v3,e1,e2,e3,ie1,cb2); 
                co3 = getco(d,it,nv,v1,v2,v3,e1,e2,e3,ie1,cb3); 

                % compute partial derivatives and store 
                partial_x_1 = decastder(d,1,b1,b2,b3,1,0,0,co1);
                partial_y_1 = decastder(d,1,b1,b2,b3,0,1,0,co1);
                partial_z_1 = decastder(d,1,b1,b2,b3,0,0,1,co1);
                GradPhi(v1(it), nq, 1) = partial_x_1;
                GradPhi(v1(it), nq, 2) = partial_y_1;
                GradPhi(v1(it), nq, 3) = partial_z_1;

                partial_x_2 = decastder(d,1,b1,b2,b3,1,0,0,co2);
                partial_y_2 = decastder(d,1,b1,b2,b3,0,1,0,co2);
                partial_z_2 = decastder(d,1,b1,b2,b3,0,0,1,co2);
                GradPhi(v2(it), nq, 1) = partial_x_2;
                GradPhi(v2(it), nq, 2) = partial_y_2;
                GradPhi(v2(it), nq, 3) = partial_z_2;

                partial_x_3 = decastder(d,1,b1,b2,b3,1,0,0,co3);
                partial_y_3 = decastder(d,1,b1,b2,b3,0,1,0,co3);
                partial_z_3 = decastder(d,1,b1,b2,b3,0,0,1,co3);
                GradPhi(v3(it), nq, 1) = partial_x_3;
                GradPhi(v3(it), nq, 2) = partial_y_3;
                GradPhi(v3(it), nq, 3) = partial_z_3;
            end 
        end 
        
        function obj = InnerProduct(obj,varargin)
            %InnerProduct compute L2 inner product matrix 
            % Usage: obj.InnerProduct()
            % Input:
            %      varargin optional arguments
            %      q: number of quadrature points for numerical
            %      integration on S2
            % Output:
            %  Basis object   
            
            p = inputParser;
            addParameter(p, 'q', 3074, @isnumeric)
            parse(p, varargin{:});
            params = p.Results;
            
            q = params.q;
            
            x = obj.Verts(:,1);
            y = obj.Verts(:,2);
            z = obj.Verts(:,3);
            
            v1 = obj.TRI(:,1);
            v2 = obj.TRI(:,2);
            v3 = obj.TRI(:,3);
            
            e1 = obj.E(:,1);
            e2 = obj.E(:,2);
            e3 = obj.E(:,3);
            
            ie1 = obj.IE(:,1);
            
            m = obj.M;
            d = obj.D;
            
            leb = getLebedevSphere(q);
            xq = leb.x;
            yq = leb.y;
            zq = leb.z;
            wq = leb.w;
            
            nv = length(v1);
            PhiQuad = zeros(m, q); 
            for nq=1:q
                % get triangle quadrature point is within 
                [it,b1,b2,b3] = sfindtri(x,y,z,v1,v2,v3,xq(nq),yq(nq),zq(nq));
                % get basis functions assocaited with triangle 
                cb1 = double(1:m == v1(it)); cb2 = double(1:m == v2(it)); cb3 = double(1:m == v3(it));
                % get local coordiantes and evaluate 
                co1 = getco(d,it,nv,v1,v2,v3,e1,e2,e3,ie1,cb1); 
                val1 = decast(d,b1,b2,b3,co1);

                co2 = getco(d,it,nv,v1,v2,v3,e1,e2,e3,ie1,cb2); 
                val2 = decast(d,b1,b2,b3,co2);

                co3 = getco(d,it,nv,v1,v2,v3,e1,e2,e3,ie1,cb3); 
                val3 = decast(d,b1,b2,b3,co3);
                % Store quadrature evaluations  
                PhiQuad(v1(it), nq) = val1;
                PhiQuad(v2(it), nq) = val2;
                PhiQuad(v3(it), nq) = val3;

            end 
            % Compute inner product matrix 
            obj.J = PhiQuad*diag(wq)*PhiQuad'; 
        end 
        
        function obj = QVInnerProduct(obj,varargin)
            %QVInnerProduct compute pair-wise inner product matrix of gradients  
            % Usage: obj.QVInnerProduct()
            % Input:
            %      varargin optional arguments
            %      q: number of quadrature points for numerical
            %      integration on S2
            % Output:
            %  Basis object  

            p = inputParser;
            addParameter(p, 'q', 3074, @isnumeric)
            parse(p, varargin{:});
            params = p.Results;
            
            q = params.q;
            
            x = obj.Verts(:,1);
            y = obj.Verts(:,2);
            z = obj.Verts(:,3);
            
            v1 = obj.TRI(:,1);
            v2 = obj.TRI(:,2);
            v3 = obj.TRI(:,3);
            
            m = obj.M;
            
            leb = getLebedevSphere(q);
            xq = leb.x;
            yq = leb.y;
            zq = leb.z;
            wq = leb.w;
            
            % Compute gradient of each basis function over all quadrature points 
            Grad_PhiQuad = obj.EvalGrad([xq, yq, zq]);
            
            % Compute \int_{S^2} [\grad\phi_i]'\grad\phi_j using quadature 
            TR = triangulation([v1, v2, v3], x, y, z); %create triangulation object for easy index
            Rmat = zeros(m,m);
            for i=1:m
                % only need to look at other basis functions in the support of \phi_i
                % get unique vertices attached to vertex i
                triRing = vertexAttachments(TR,i);
                vertRing = TR.ConnectivityList(triRing{:},:);
                vertRing_unique = unique(vertRing)';

                for j=vertRing_unique
                    grad_phi_iq = squeeze(Grad_PhiQuad(i,:,:));
                    grad_phi_jq = squeeze(Grad_PhiQuad(j,:,:));
                    for nq=1:q
                        Rmat(i,j) = Rmat(i,j) + wq(nq)*grad_phi_iq(nq,:)*grad_phi_jq(nq,:)';
                    end 
                end 
            end 
            obj.R = Rmat;
        end 
        
        function plot(obj)
            % plot plot the tesselation
            %
            % Usage: obj.plot()
            x = obj.Verts(:,1);
            y = obj.Verts(:,2);
            z = obj.Verts(:,3);

            ie1 = obj.IE(:,1);
            ie2 = obj.IE(:,2);

            srendtri(x,y,z,ie1,ie2)
        end 
        
    end
end

