function [mat, supp] = heinestieltjes(Q, r, n)
    % r = Fuchs index
    % Q = k poly in cell array with real coefficients (deg 0 to k + r)
    % n = degree of the polynomial solution
    
    A = zeros(n + 1, n + r + 1);
    for i = 0:n
        for j = 0:n + r
            A(n + 1 - i, n + r + 1 - j) = L(i, j, Q);
        end
    end

    mat = cell(r + 2, 1);
    mat{1} = A.';
    for i = 0:r
        mat{i + 2} = [zeros(n + 1, i) eye(n + 1) zeros(n + 1, r - i)].';
    end
    supp = [zeros(1, r + 1); eye(r + 1)];
end

function l = L(p, q, Q)
    k = length(Q);
    l = 0;
    for r = 1:k
        l = l + (factorial(p)/factorial(p - r))*Q{r}(q - p + r + 1);
    end
end
