function [mat, supp] = heinestieltjes(Q, n)
    % Q = k poly in cell array with real coefficients (deg 0 to k + r)
    % n = degree of the polynomial solution

       
    % r = Fuchs index: now detected from polynomials

    k = length(Q);
    r = -k;
    for i = 1:k
        degQi = length(Q{i}) - 1;
        r = max(r, degQi - i);
    end
    r
    
    A = zeros(n + 1, n + r + 1);
    for i = 0:n
        for j = 0:i + r
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
    for i = 1:k
        if p >= i
            degQi = length(Q{i}) - 1;
            coef = q - p + i;
            if degQi >= coef
                if coef > 0
                    l = l + (factorial(p)/factorial(p - i))*Q{i}(coef + 1);
                end
            end
        end
    end
end
