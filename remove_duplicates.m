function B = remove_duplicates(A,tol)

% REMOVE_DUPLICATES removes multiple rows from matrix A
% 
% B = remove_duplicates(A,tol) returns just unique of rows of A
% using relative tolerance tol

% Bor Plestenjak 2026

if nargin<2
    tol = 1e-6;
end

n = size(A,1);
m = size(A,2);
for j = 1:n
    elem = A(j,:);
    if min(vecnorm(double(A([1:j-1 j+1:end],:) - elem),2,2))/(1+norm(elem)) < tol
        A(j,:) = Inf*ones(1,m);
    end
end

pos = find(vecnorm(double(A),2,2)<Inf);
B = A(pos,:);
