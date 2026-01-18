function [B,lambda,X] = initial_rmep3(n)

%INITIAL_RMEP3   Initial linear 3-parameter rectangular MEP for the homotopy
% [B,lambda,X] = initial_rmep3(n) returns matrices, eigenvalues 
% and eigenvectors of a RMEP
%
% (B0 + lambda(1) B1 + lambda(2) B2 + lambda(3)B3) x = 0 
%
% Input:
%   - n: number of columns in matrices B
%
% Output:
%   - B : cell {B0,B1,B2,B3}
%   - lambda : matrix m x k, each row is an eigenvalue
%   - X : matrix n x m with right eigenvectors

% Bor Plestenjak, 2026

a = randn(n,1)+1i*randn(n,1);
b = randn(n,1)+1i*randn(n,1);

B0 = [diag(a);zeros(2,n)] + [zeros(2,n); diag(b)];
B1 = [eye(n);zeros(2,n)];
B2 = [zeros(1,n); eye(n); zeros(1,n)];
B3 = [zeros(2,n); eye(n)];

B = {B0,B1,B2,B3};

nsol = n*(n+1)*(n+2)/6;
lambda = zeros(nsol,3);
X = zeros(n,nsol);
ind_col = 0;
z0 = zeros(n,1);
for i = 1:n
    ind_col = ind_col + 1;
    z = z0;
    z(i) = 1;
    lambda(ind_col, :)= [-a(i) 0 -b(i)];
    X(:,ind_col) = z;
    for j = i+1:n
        lm = -a(i);
        um = -b(j);
        A = B0(i+1:j+1,i:j) + lm*B1(i+1:j+1,i:j) + um*B3(i+1:j+1,i:j);
        [Z,D] = eig(A);
        d = diag(D);
        for r = 1:length(d)
            ind_col = ind_col + 1;
            z = z0;
            z(i:j) = Z(:,r);
            lambda(ind_col, :)= [lm -d(r) um];
            X(:,ind_col) = z/norm(z);
        end
    end
end