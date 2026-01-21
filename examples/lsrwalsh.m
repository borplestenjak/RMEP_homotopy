function [mat,supp] = lsrwalsh(y,n)
    %LSR   Least-squares realization problem via Walsh' theorem.
    %   mat = LSRWALSH(y,n) creates the coefficient matrices of the
    %   multiparameter eigenvalue problem that solves the n-th order least-
    %   squares realization problem for the output data in the vector y.
    %
    %   [mat,supp] = LSRWALSH(...) also gives the support of the 
    %   coefficient matrices.
        
    % Copyright (c) 2024 - Christof Vermeersch
    %
    % This function implements the construction of the numerical
    % coefficient matrices as described in [1].
    %
    % [1] Lagauw, S., Vanpoucke, L., and De Moor, B. (2024) Exact
    % Characterizationof the Global Optima of Least Squares Realization of
    % Autonomous LTI Models as a Multiparameter Eigenvalue Problem, Proc. 
    % of the 22nd European Control Conference (ECC), Stockholm, Sweden.

    % Define the (symbolic) variables:
    N = length(y);
    a = sym("a",[n 1]);
    for i = 1:n
        assume(a(i),"real")
    end

    % Create the auxiliary matrices:
    t = [fliplr(a') 1 zeros(1,N-n-1)];
    Ta = sym('Ta',[N-n,N]);
    for i = 1:N-n
        Ta(i,:) = circshift(t,i-1);
    end
    Tb = Ta(1:N-2*n,1:N-n);

    % Construct matrix polynomial:
    A = [Ta*y Ta*Ta'*Tb'];

    % Derive the numerical coefficient matrices:
    [monomials,supp] = powers(a);
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
    %   All monomials of total degree up to three are required by the main
    %   function.

    n = length(a);
    I = fliplr(eye(n,n));
    supp = zeros(nchoosek(n+3,3),n);
    supp(2:n+1,:) = I;
    endPoint = 1;
    blockLength = ones(n,1);
    blockPos = [0; cumsum(blockLength)];
    for di = 2:3
        blockLength = cumsum(blockLength,"reverse");
        blockPos = [0; cumsum(blockLength)];
        endPoint = endPoint + blockLength(1);
        for ni = n:-1:1
            supp(endPoint+blockPos(ni)+1:endPoint+blockPos(ni+1),:) = ...
                supp(endPoint-blockLength(ni)+1:endPoint,:) + I(ni,:);
        end
    end
    monomials = prod(a'.^supp,2);
end
