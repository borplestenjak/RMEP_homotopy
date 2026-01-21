function [mat,supp] = arma(y,na,nc)
    %ARMA   Autoregressive moving-average problem.
    %   mat = ARMA(y,na,nc) creates the coefficient matrices of the
    %   multiparameter eigenvalue problem that solves the autoregressive
    %   moving-average problem with orders na and nc for the output data 
    %   in the vector y.
    %
    %   [mat,supp] = ARMA(...) also gives the support of the coefficient
    %   matrices.
        
    % Copyright (c) 2024 - Christof Vermeersch
    %
    % This function implements the construction of the numerical
    % coefficient matrices as described in [1].
    %
    % [1] Vermeersch, C. and De Moor, B. (2019) Globally Optimal
    % Least-Squares ARMA Model Identification Is an Eigenvalue Problem,
    % IEEE Control Systems Letters, vol. 3, no. 4, p. 1063-1067.

    % Define the (symbolic) variables:
    N = length(y);
    M = N-na;
    n = na+nc;
    a = sym("a",[na 1]);
    for i = 1:na
        assume(a(i),"real")
    end
    c = sym("c",[nc 1]);
    for i = 1:nc
        assume(c(i),"real")
    end

    % Create the auxiliary matrices:
    t = [fliplr(a') 1 zeros(1,M-1)];
    Ta = sym('Ta',[M,N]);
    for i = 1:M
        Ta(i,:) = circshift(t,i-1);
    end
    t = [fliplr(c') 1 zeros(1,M-1)];
    Tc = sym('Ta',[M,M+nc]);
    for i = 1:M
        Tc(i,:) = circshift(t,i-1);
    end
    Dc = Tc*Tc';
    
    % Construct matrix polynomial:
    A = sym("A",[M*(n+1)+n M*(n+1)+1]);
    A(:,:) = zeros(M*(n+1)+n,M*(n+1)+1);
    for i = 1:na
        A(i,(i-1)*M+1:i*M) = y.'*Ta.';
        A(i,M*n+1:M*(n+1)) = (diff(Ta,a(i))*y).';
    end
    for i = 1:nc
        A(i+na,(na+i-1)*M+1:(na+i)*M) = y.'*Ta.';
    end
    for i = 1:na
        A(n+(i-1)*M+1:n+i*M,(i-1)*M+1:i*M) = Dc;
        A(n+(i-1)*M+1:n+i*M,end) = diff(Ta,a(i))*y;
    end
    for i = 1:nc
        A(n+(na+i-1)*M+1:n+(na+i)*M,(na+i-1)*M+1:(na+i)*M) = Dc;
        A(n+(na+i-1)*M+1:n+(na+i)*M,M*n+1:M*(n+1)) = diff(Dc,c(i));
    end
    A(end-N+na+1:end,M*n+1:M*(n+1)) = Dc;
    A(end-N+na+1:end,end) = Ta*y;

    % Derive the numerical coefficient matrices:
    [monomials,supp] = powers(a,c);
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
function [monomials,supp] = powers(a,c)
    %POWERS   Monomials and support.
    %   [monomials,supp] = POWERS(a,c) returns the symbolic monomials and 
    %   the (numerical) support for the different powers of the variables 
    %   in a and c. All monomials in c of total degree up to two, but only 
    %   the linear monomials in a, are required by the main function.

    n = length(c);
    I = fliplr(eye(n,n));
    supp = zeros(nchoosek(n+2,2),n);
    supp(2:n+1,:) = I;
    blockLength = cumsum(ones(n,1),"reverse");
    blockPos = [0; cumsum(blockLength)];
    for ni = n:-1:1
        supp(n+blockPos(ni)+2:n+blockPos(ni+1)+1,:) = ...
            supp(n-blockLength(ni)+2:n+1,:) + I(ni,:);
    end
    monomials = [prod(c'.^supp,2); a];
    supp = [zeros(size(supp,1),length(a)) supp; ...
        eye(length(a)) zeros(length(a),length(c))];
end