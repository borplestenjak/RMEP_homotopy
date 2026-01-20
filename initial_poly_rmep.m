function [B,lambda,X] = initial_poly_rmep(n,k,deg)

%INITIAL_POLY_RMEP   Initial polynomial rectangular MEP for the homotopy
% [B,lambda,X] = initial_rmep(n,k) returns matrices, eigenvalues 
% and eigenvectors of a RMEP
%
% (B0 + lambda(1)^deg(1) B1 + ... + lambda(k)^deg(k) Bk) x = 0 
%
% Input:
%   - n: number of columns in matrices B
%   - k: number of parameters
%   - deg: degrees of monomials
%
% Output:
%   - B : cell {B0,B1,...,Bk}
%   - lambda : matrix m x k, each row is an eigenvalue
%   - X : matrix n x m with right eigenvectors

% Bor Plestenjak, 2026

% we construct linear RMEP first
[B,mu,Xlin] = initial_rmep(n,k);
nLin = size(mu,1);

% and convert its eigenvalues by computing alls possible roots
lambda = [];
X = [];

for i = 1:nLin
    M = [];
    for j = 1:k
        fi = (0:deg(j)-1)*2*pi/deg(j);
        omega = cos(fi) + 1i*sin(fi);
        lam = omega.'*mu(i,j)^(1/deg(j));
        if j == 1
            M = lam;
        else
            M = [kron(M,ones(deg(j),1)) kron(ones(size(M,1),1),lam)];
        end
    end
    lambda = [lambda; M];
    X = [X kron(Xlin(:,i),ones(1,size(M,1)))];
end


