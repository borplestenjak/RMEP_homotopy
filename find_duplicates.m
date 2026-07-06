function ind = find_duplicates(M,tol)

% FIND_DUPLICATES returns indices of equal rows in matrix A
% 
% ind = find_duplicates(A,tol) returns indices of rows in A that are not 
% unique using relative tolerance tol

% Bor Plestenjak 2026

    if nargin<2
        tol = 1e-6;
    end
    
    % we use sketching with a random direction to filter candidates
    g = randn(size(M,2),1);
    g = g/norm(g);
    z = M*g;
    
    cand = find_duplicates_complex_vec(z,tol);
    
    A = M(cand,:);
    
    n = size(A,1);
    isdup = false(n,1);
    
    scale = 1 + vecnorm(double(A),2,2);
    
    for j = 1:n
        d = vecnorm(double(A - A(j,:)),2,2);
        d(j) = inf;
    
        if min(d) < tol*scale(j)
            isdup(j) = true;
        end
    end
    
    indA = find(isdup).';
    ind = cand(indA);

end

% helper function finds indices with close elements in a complex vector

function ind = find_duplicates_complex_vec(z,tol)

    n = numel(z);
    
    [~,perm] = sort(real(z));
    zs = z(perm);
    
    dup = false(n,1);
    
    for j = 1:n
        k = j+1;
        thresh = tol*(1+abs(zs(j)));
    
        while k <= n && abs(real(zs(k))-real(zs(j))) < thresh
            if abs(zs(k)-zs(j)) < thresh
                dup(j) = true;
                dup(k) = true;
            end
            k = k+1;
        end
    end
    
    ind = sort(perm(dup)).';

end