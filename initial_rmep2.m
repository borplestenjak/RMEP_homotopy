function [B,lambda,X] = initial_rmep2(n)

%INITIAL_RMEP2   Initial linear 2-parameter rectangular MEP for the homotopy
% [B,lambda,X] = initial_rmep2(n) returns matrices, eigenvalues 
% and eigenvectors of a RMEP
%
% (B0 + lambda(1) B1 + lambda(2) B2) x = 0 
%
% Input:
%   - n: number of columns in matrices B
%
% Output:
%   - B : cell {B0,B1,B2}
%   - lambda : matrix m x k, each row is an eigenvalue
%   - X : matrix n x m with right eigenvectors

% Bor Plestenjak, 2026

a = randn(n,1)+1i*randn(n,1);
b = randn(n,1)+1i*randn(n,1);

B0 = [diag(a);zeros(1,n)] + [zeros(1,n); diag(b)];
B1 = [eye(n);zeros(1,n)];
B2 = [zeros(1,n); eye(n)];

B = {B0,B1,B2};

nsol = n*(n+1)/2;
lambda = zeros(nsol,2);
X = zeros(n,nsol);
ind_col = 0;
z0 = zeros(n,1);
for i = 1:n
    for j = i:n
        lm = -a(i);
        um = -b(j);
        ind_col = ind_col + 1;
        lambda(ind_col, :)= [lm um];
        z = z0;
        z(i) = 1;
        for r = i+1:j
            z(r) = -(um + b(r-1))*z(r-1)/(a(r)+lm);
        end
        X(:,ind_col) = z/norm(z);
    end
end