function [Vtrajectory] = power_iterations(M, v0, maxiter, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function to  leading eigenvector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inputs:
%M - the m x m matrix  
%v0 - mx1 initial guess 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%outputs:
%v - m x 1: leading eigenvectors 

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
    Vtrajectory(:,i+1) = M*Vtrajectory(:,i)/norm(M*Vtrajectory(:,i));
    delta_residual_norm = norm(Vtrajectory(:,i+1) - Vtrajectory(:,i));
    i = i + 1;
end 
Vtrajectory = Vtrajectory(:,1:i);