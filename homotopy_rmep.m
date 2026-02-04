function [lambda,X,tcell,ycell,stat] = homotopy_rmep(A,B,Lambda0,X0,opts)

% HOMOTOPY_RMEP   Solve a linear rectangular multiparameter eigenvalue 
% problem using homotopy
%
% [lambda,X,tn,yn,stat] = homotopy_rmep(A,B,Lambda0,X0,opts) returns 
% homogeneous eigenvalues and eigenvectors of a rectangular multiparameter 
% eigenvalue problem (RMEP)
%
% (A{1} + lambda(1) A{2} x + ... + lambda(k) A{k+1}) x = 0 
%
% using initial RMEP
%
% (B{1} + lambda(1) B{2} x + ... + lambda(k) B{k+1}) x = 0 
%
% with known eigenvalues Lambda0 and eigenvectors X0
%
% Input:
%   - A,B : cell arrays of size k+1, each cell is an (n+k-1) x n matrix 
%   - Lambda0: rows are eigenvalues of the RMEP B (matrix of size m x (k+1)),
%         where m = nchoosek(n+k-1,k)
%   - X0: columns are eigenvectors of the RMEP B (matrix of size n x m)
%   - opts : options (for the list of options see homotopy_poly_rmep)
%
% Output:
%   - lambda : matrix m x (k+1), each row is a homogeneous eigenvalue
%   - X : matrix n x m with right eigenvectors
%   - tcell: cell with all time steps 
%   - ycell: cell with all intermediate steps, first n elements represent eigenvector,
%         last k+1 elements represent eigenvalue
%   - stat: several statistics

% Bor Plestenjak 2026

% Validate number of input parameters
narginchk(4, 5);
if nargin < 5, opts = []; end

[lambda,X,tcell,ycell,stat] = homotopy_poly_rmep(A,[],B,[],Lambda0,X0,opts);