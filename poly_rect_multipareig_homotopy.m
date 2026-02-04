function [lambda,X,lambdaT,XT] = poly_rect_multipareig_homotopy(A,suppA,opts)

%POLY_RECT_MULTIPAREIG_HOMOTOPY  Solve a polynomial rectangular 
% multiparameter eigenvalue problem using homotopy
%
% [lambda,X,lambdaT,XT] = poly_rect_multipareig_homotopy(A,suppA,degB,opts) 
% returns eigenvalues and eigenvectors of a polynomial rectangular 
% multiparameter eigenvalue problem (RMEP)
%
% A(lambda) x = sum(j=1)^r1 (lambda.^suppA_r(j,:)*Aj) x = 0 
% 
% Input:
%   - A: cell arrays containing coefficient matrices of a polynomial RMEP, 
%         each cell is an (n+k-1) x n matrix
%   - suppA: matrices with rows that contain degrees of monomials coefficient 
%         if suppA is empty, we assume that RMEP is linear
%   - opts: options
%
% Options in opts:
%   - filter_res (1e-8): return only eigenpairs with small residual
%       (use when homogeneous RMEPS has positive-dimensional solution, or
%        in case of unwanted inifinite eignvalues), set to zero to return
%        all eigenpairs
%   - degB: degrees used for the initial problem B
%   - ali options available in homotopy_poly_rmep
%
% Output:
%   - lambda : matrix m1 x k, each row is an eigenvalue
%   - X : matrix n x m1 with right eigenvectors for lambda
%   - lambdaT : matrix m2 x (k+1), each row is a homogenous eigenvalue
%     (might include additional eigenvalues if filter_res>0) 
%   - XT : matrix n x m2 with right eigenvectors for lambdaT
      
% Bor Plestenjak 2026

if nargin < 2, suppA = []; end
if nargin < 3, opts = []; end

class_t = superiorfloat(A{:});

if isfield(opts,'filter_res'), filter_res = opts.filter_res;  else,  filter_res = sqrt(eps(class_t));   end
if isfield(opts,'degB'),       degB = opts.degB;              else,  degB = [];                         end

if isempty(suppA)
    k = numel(A) - 1;
else
    k = size(suppA,2);   % number of parameters
end

n = size(A{1},2);

% for each of the variables we check maximal degree of the monomials that
% include the variable
if isempty(degB)
    degrows = sum(suppA,2);
    if isempty(suppA)
        degB = ones(1,k);
    else
        for j = 1:k 
            indrow = find(suppA(:,j)>0);
            degB(1,j) = max(degrows(indrow));
        end
    end
end

[B,Lambda0,X0] = initial_poly_rmep(n,k,degB);
suppB = [zeros(1,k); diag(degB)]; 
[lambdaT, XT] = homotopy_poly_rmep(A,suppA,B,suppB,Lambda0,X0,opts);
% homotopy is solving homogeneous problem, solutions might include infinite eigenvalues

% conversion back from homogeneous eigenvalues
lambdaN = lambdaT(:,2:end)./lambdaT(:,1);

if filter_res>0
    res = [];
    for j = 1:length(lambdaN)
        W = eval_rmep(A,suppA,lambdaN(j,:));
        res(j,:) = norm(W*XT(:,j));
    end
    ind = find(res<filter_res);
else
    ind = 1:length(lambdaN);
end

lambda = lambdaN(ind,:);
X = XT(:,ind);