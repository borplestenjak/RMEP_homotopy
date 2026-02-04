function ans = is_in_set(A,x,tol)

% IS_IN_SET check if row is in a given matrix
% 
% ans = is_in_set(A,x,tol) checks (using tolerance tol) if x is close
% enough to a row of matrix A
%
% Bor Plestenjak 2026

class_t = superiorfloat(A,x);

if nargin<3
    tol = 1e-6;
end

m = size(A,1);
if m==0
    ans = 0;
else
    if strcmp(class_t,'mp')
        ans = min(vecnorm(double(A-x),2,2)) <= tol*(1+norm(x));
    else
        ans = min(vecnorm(A-x,2,2)) <= tol*(1+norm(x));
    end
end
