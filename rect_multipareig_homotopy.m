function [lambda,X,lambdaT,XT,stat] = rect_multipareig_homotopy(A,opts)

%RECT_MULTIPAREIG_HOMOTOPY  Solve a linear rectangular multiparameter 
% eigenvalue problem using homotopy
%
% [lambda,X,lambdaT,XT] = rect_multipareig_homotopy(A,opts) returns 
% eigenvalues and eigenvectors of a linear RMEP
%
% A{1} x + lambda(1) A{2} x + ... + lambda(k) A{k+1} x = 0 
% 
% Input:
%   - A : cell array of size k+1 of matrices A{i}, all matrices have to be 
%         rectangular matrices of the same size (n+k-1) x n
%   - opts: options
%
% Options in opts:
%   - filter_res (1e-8): return only eigenpairs with small residual
%       (use when homogeneous RMEPS has positive-dimensional solution, or
%        in case of unwanted inifinite eignvalues), set to zero to return
%        all eigenpairs
%   - ali options available in homotopy_poly_rmep
%
% Output:
%   - lambda : matrix m1 x k, each row is an eigenvalue
%   - X : matrix n x m1 with right eigenvectors
%   - lambdaT : matrix m2 x (k+1), each row is a homogenous eigenvalue
%     (might include additional eigenvalues if filter_res>0) 
%   - XT : matrix n x m2 with right eigenvectors for lambdaT
      
% Bor Plestenjak 2026

if nargin < 2, opts = []; end

class_t = superiorfloat(A{:});

if isfield(opts,'filter_res'),      filter_res = opts.filter_res;           else,  filter_res = sqrt(eps(class_t));            end

k = numel(A) - 1;
n = size(A{1},2);

[B,Lambda0,X0] = initial_rmep(n,k);
[lambdaT, XT, ~, ~, stat] = homotopy_rmep(A,B,Lambda0,X0,opts);
% homotopy is solving homogeneous problem, solutions might include infinite eigenvalues

% conversion back from homogeneous eigenvalues
lambdaN = lambdaT(:,2:end)./lambdaT(:,1);
neig = size(lambdaN,1);
suppA = [zeros(1,k); eye(k)];

if filter_res>0
    res = [];
    for j = 1:neig
        W = eval_rmep(A,suppA,lambdaN(j,:));
        res(j,:) = norm(W*XT(:,j));
    end
    ind = find(res<filter_res);
else
    ind = 1:length(neig);
end

lambda = lambdaN(ind,:);
X = XT(:,ind);

