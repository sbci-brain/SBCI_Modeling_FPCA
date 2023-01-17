function [Vtrajectory] = generalized_power(M, J, v0, eta, maxiter, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to compute sparse leading eigenvector
%algorithm proposed in: DOI:10.5555/1756006.1756021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inputs:
%M - the m x m matrix  
%v0 - mx1 initial guess 
%eta - number of non-zero components 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%outputs:
%v - m x 1: leading sparse eigenvectors 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%symmetry check 
if ~issymmetric(M)
    M = (M + M')/2;
end 

i = 1; delta_residual_norm = inf;
Vtrajectory = zeros(size(M,1), maxiter);
Vtrajectory(:,1) = v0;
while (i < maxiter) && (delta_residual_norm > tol)
    Vtrajectory(:,i+1) = M*Vtrajectory(:,i);
    [out, ix_i] = sort(abs(Vtrajectory(:,i+1)), 'descend');
    Vtrajectory(ix_i(eta+1:end),i+1) = 0;
    Jnorm_i = sqrt(Vtrajectory(:,i+1)'*J*Vtrajectory(:,i+1));
    Vtrajectory(:,i+1) = Vtrajectory(:,i+1)/Jnorm_i;
    delta_residual_norm = norm(Vtrajectory(:,i+1) - Vtrajectory(:,i));
    i = i + 1;
end 
Vtrajectory = Vtrajectory(:,1:i);