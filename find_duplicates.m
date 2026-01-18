function ind = find_duplicates(A,tol)

% FIND_DUPLICATES returns indices of equal rows in matrix A
% 
% ind = find_duplicates(A,tol) returns indices of rows in A that are not 
% unique using relative tolerance tol

% Bor Plestenjak 2026

if nargin<2
    tol = 1e-6;
end

ind = [];

n = size(A,1);
m = size(A,2);
for j = 1:n
    elem = A(j,:);
    if min(vecnorm(A([1:j-1 j+1:end],:) - elem,2,2))/(1+norm(elem)) < tol
        ind = [ind j];
    end
end
