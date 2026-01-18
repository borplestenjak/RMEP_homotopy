function [B,lambda,X] = initial_rmep(n,k)

%INITIAL_RMEP   Initial linear rectangular MEP for the homotopy
% [B,lambda,X] = initial_rmep(n,k) returns matrices, eigenvalues 
% and eigenvectors of a RMEP
%
% (B0 + lambda(1) B1 + ... + lambda(k) Bk) x = 0 
%
% Input:
%   - n: number of columns in matrices B
%   - k: number of parameters
%
% Output:
%   - B : cell {B0,B1,...,Bk}
%   - lambda : matrix m x k, each row is an eigenvalue
%   - X : matrix n x m with right eigenvectors

% Bor Plestenjak, 2026

if k==2
    % fast construction for 2-parameter problems
    [B,lambda,X] = initial_rmep2(n);
elseif k==3
    % fast construction for 3-parameter problems
    [B,lambda,X] = initial_rmep3(n);
else
    % recursive construction for problems with more than 3 parameters
    [B,lambda,X] = initial_rmep_rec(n,k);
end

