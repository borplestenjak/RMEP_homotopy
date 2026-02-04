function [gm,s,W,B] = condeig_rmep(A,suppA,lambda,affine,tol)

%CONDEIG_RMEP Geometric multiplicity and condition number of an eigenvalue 
% of a polynomial rectangular multiparameter eigenvalue problem 
%
% A(lambda) x = sum(j=1)^r1 (lambda.^suppA_r(j,:)*Aj) x = 0 
%
% If suppA is empty, we assume that A is a linear RMEP
%
% usage: [gm,s,W,B] = condeig_rmep(A,suppA,lambda,affine,tol) 
%
% Input:
%   - A: cell array containing coefficient matrices of a polynomial RMEP, 
%         each cell is an (n+k-1) x n matrix
%   - suppA: matrices with rows that contain degrees of monomials coefficient 
%         if suppA is empty, we assume that the RMEP is linear
%   - lambda: eigenvalue (row of length k)
%   - affine: is problem in affine form (true:default) or homogeneous (false)
%   - tol: tolerance for multiple eigenvalue 
%
% Output:
%   - gm: geometric multiplicity of eigenvalue lambda (0 if not an eigenvalue)
%   - s : condition number of the eignvalue (=cond(B))
%   - W : matrix A(lambda)
%   - B : for a geometrically simple eigenvalue, this is a k x k matrix, such
%         that B(i,j) = y_i'*part(A,lambda_j)*x, where y_i is a vector from
%         orthonormal basis of ker(A(lambda)'), part(A,lambda_j) is a partial 
%         derivative of A(lambda) computed at lambda, and x is an eigenvector

% Bor Plestenjak 2026

class_t = superiorfloat(A{:});

if nargin<4
    affine = true;
end

[m,n] = size(A{1});
n_cellA = numel(A);
if isempty(suppA)
    if affine
        k = n_cellA - 1;
        suppA = [zeros(1,k); eye(k)];
    else
        k = n_cellA;
        suppA = eye(k);
    end
else
    k = size(suppA,2);
end

W = eval_rmep(A,suppA,lambda);
if nargin<5
    tol = sqrt(eps)*norm(W);
end
s = svd(W);
gm = n - rank(W,tol);
if gm == 1
    [U,S,V] = svd(W);
    x = V(:,end);
    if affine
        Y = U(:,end-k+1:end);
    else
        Y = U(:,end-k+2:end);
    end
    Ca = zeros(m,n_cellA,class_t);  
    for j = 1:n_cellA
        Ca(:,j) = A{j}*x;
    end
    DW = eval_rmep_der(Ca,suppA,lambda);
    B = Y'*DW;
    s = cond(B);
else
    B = [];
    s = Inf;
end
