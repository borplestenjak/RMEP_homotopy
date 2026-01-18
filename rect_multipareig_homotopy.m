function [lambda,X] = rect_multipareig_homotopy(A,opts)

%RECT_MULTIPAREIG_HOMOTOPY  Solve a linear rectangular multiparameter 
% eigenvalue problem using homotopy
%
% [lambda,X] = rect_multipareig_homotopy(A,opts) returns eigenvalues and 
% eigenvectors of a rectangular multiparameter eigenvalue problem
%
% A{1} x + lambda(1) A{2} x + ... + lambda(k) A{k+1} x = 0 
% 
% Input:
%   - A : cell array of size k+1 of matrices A{i}, all matrices have to be 
%         rectangular matrices of the same size (n+k-1) x n
%   - opts : options 
%
% Options in opts: (for more see homotopy_rmep)
%   - dist_tol (1e-6): relative tolerance for checking distinct solutions
      
% Bor Plestenjak 2026

if nargin < 2, opts = []; end

if isfield(opts,'dist_tol'),  dist_tol = opts.dist_tol;  else,  dist_tol = 1e-6;  end

k = size(A,2) - 1;
n = size(A{1},2);

run = 1;    
[B,Lambda0,X0] = initial_rmep(n,k);
[lambdaT, X, tn, yn, stat] = homotopy_rmep(A,B,Lambda0,X0,opts);
lambda = lambdaT(:,2:end)./lambdaT(:,1);
ind = find_duplicates(lambda,dist_tol);

while length(ind)>0
    run = run + 1;
    fprintf('Another attempt %d\n',run)
    [B,Lambda0,X0] = initial_rmep(n,k);
    [lambdaT, X, tn, yn] = homotopy_rmep(A,B,Lambda0,X0,opts);
    lambda = lambdaT(:,2:end)./lambdaT(:,1);
    ind = find_duplicates(lambda,dist_tol);
end    


