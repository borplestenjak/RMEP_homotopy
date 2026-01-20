function Y = eval_rmep(A,supp,x,Y)
    
    % EVAL_RMEP evaluate the polynomial multivariate matrix in the 
    % standard monomial basis.
    %
    %   Y = EVAL_RMEP(A,supp,x,Y) evaluates the rectangular polynomial 
    %   matrix given by its representation in the standard monomial basis.
    %
    %   Input argument:
    %       A: cell of size r with coefficient matrices
    %       supp: matrix r x k with monomial degrees 
    %       x: point (size 1 x k) to evaluate
    %       Y: (optional) zero matrix of the same same as matrices in A 
    %
    %   Output argument:
    %       Y = value of A(x)
    
    % Based on evalmon from MacaulayLab by Christof Vermeersch
    
    % Bor Plestenjak 2026
    
    if nargin<4
        % When called in homotopy_poly_rmep it is faster to provide
        % zero matrix as argument
        [m,n] = size(A{1});
        Y = zeros(m,n);
    end

    x = x(:).'; % must be a row
    coef = prod(x.^supp,2);
    for j = 1:size(supp,1)
        Y = Y + coef(j)*A{j};
    end
