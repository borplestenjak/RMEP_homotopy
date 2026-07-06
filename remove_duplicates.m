function B = remove_duplicates(M,tol)

% REMOVE_DUPLICATES removes multiple rows from matrix A
% 
% B = remove_duplicates(A,tol) returns just unique of rows of A
% using relative tolerance tol

% Bor Plestenjak 2026

if nargin<2
    tol = 1e-6;
end

n = size(M,1);
ind_mult = find_duplicates(M,tol);
indices = ones(1,n);
indices(ind_mult) = 0;
pos = find(indices==1);
B = M(pos,:);

% A contains only multiple rows, we select just unique rows
A = M(ind_mult,:);
n = size(A,1);
m = size(A,2);
for j = 1:n
    elem = A(j,:);
    if min(vecnorm(double(A([1:j-1 j+1:end],:) - elem),2,2))/(1+norm(elem)) < tol
        A(j,:) = Inf*ones(1,m);
    end
end

pos = find(vecnorm(double(A),2,2)<Inf);
B = [B; A(pos,:)];
