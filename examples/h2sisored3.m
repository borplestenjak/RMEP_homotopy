function [mat,supp] = h2sisored3(an,bn,r)
    %H2SISORED3   H2-norm model order reduction problem.
    %   mat = H2SISORED3(an,bn,r) creates the coefficient matrices of the
    %   multiparameter eigenvalue problem that solves the H2-norm model 
    %   order reduction problem of an n-th order model to an r-th order 
    %   model. It uses the so-called "third (Walsh) approach" from [1]. The 
    %   n-th order model is given by its transfer function b(1)*s^{n-1} + 
    %   b(1)*s^{n-2} + ... + b(n) / s^n + a(1)*s^{n-1} + ... + a(n).
    %
    %   [mat,supp] = H2SISORED3(...) also gives the support of the 
    %   coefficient matrices.
        
    % Copyright (c) 2024 - Christof Vermeersch
    %
    % This function implements the construction of the numerical
    % coefficient matrices as described in [1].
    %
    % [1] Lagauw, S., Agudelo, O. M., and De Moor, B. (2023) Globally
    % Optimal SISO H2-norm Model Reduction Using Walsh's Theorem, IEEE
    % Control System Letters, vol. 7, p. 1970-1675.

    % Define the (symbolic) variables:
    n = length(an);
    ar = sym("ar",[r 1]);
    for i = 1:r
        assume(ar(i),"real")
    end

    % Create the auxiliary matrices:
    t = [1 ar' zeros(1,n-1)];
    Ta = sym('Ta',[n,n+r]);
    for i = 1:n
        Ta(i,:) = circshift(t,i-1);
    end
    t = [1 an' zeros(1,r-1)];
    Tb = zeros(r,n+r);
    for i = 1:r
        Tb(i,:) = circshift(t,i-1);
    end
    t = [1 ar' zeros(1,n-1)];
    Tc = sym('Tc',[n,n+r]);
    for i = 1:n
        Tc(i,:) = circshift(t,i-1);
    end
    Td = Tc(1:n-r,1:n);

    % Construct matrix polynomial:
    A = [Ta.'*bn Tb.' Tc.'*Td.'];

    % Derive the numerical coefficient matrices:
    [monomials,supp] = powers(ar);
    mat = cell(length(monomials),1);
    for i = 1:length(monomials)
        mat{i} = zeros(size(A));
    end
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            [coefs,mons] = coeffs(A(i,j));
            for k = 1:length(mons)
                idx = find(mons(k) == monomials);
                temp = mat{idx};
                temp(i,j) = coefs(k);
                mat{idx} = temp;
            end
        end
    end
end

% Define the required subfunction:
function [monomials,supp] = powers(a)
    %POWERS   Monomials and support.
    %   [monomials,supp] = POWERS(a) returns the symbolic monomials and the
    %   (numerical) support for the different powers of the variables in a.
    %   All monomials of total degree up to two are required by the main
    %   function.

    n = length(a);
    I = fliplr(eye(n,n));
    supp = zeros(nchoosek(n+2,2),n);
    supp(2:n+1,:) = I;
    blockLength = n:-1:1;
    blockPos = [0; cumsum(blockLength')];
    for i = n:-1:1
        supp(n+blockPos(i)+2:n+blockPos(i+1)+1,:) = ...
            supp(n-blockLength(i)+2:n+1,:) + I(i,:);
    end
    monomials = prod(a'.^supp,2);
end