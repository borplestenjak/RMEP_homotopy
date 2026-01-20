function Y = eval_rmep_der(A,supp,x,coef)
    
    % EVAL_RMEP_DER evaluates products of all partial derivatives 
    % of a multivariate column A(lambda) in the standard monomial basis.
    %
    %   Y = EVAL_RMEP_DER(A,supp,x) evaluates the matrix Y whose columns
    %   are partial derivatives of A(lambda) 
    %
    %   Input argument:
    %       A: matrix of size m x r with coefficient columns
    %       supp: matrix r x k with monomial degrees 
    %       x: point (size 1 x k) to evaluate
    %       coef: (optional) supply empty matrix of size r x k
    %   Output argument:
    %       Y = matrix m x k with partial derivatives

    % Based on evalmon from MacaulayLab by Christof Vermeersch
    
    % Bor Plestenjak 2026
    
    [r, k] = size(supp);
    if nargin < 4
        % When called in homotopy_poly_rmep it is faster to provide
        % zero matrix as argument
        coef = zeros(r,k);
    end
    
    x = x(:).'; % must be a row

    for j=1:k
        der_supp = supp;
        fac = der_supp(:,j);
        der_supp(:,j) = max(der_supp(:,j)-1,0); 
        coef(:,j) = fac.*prod(x.^der_supp,2);  
    end
    Y = A*coef;
