function [value] = evalmon_der(coef,supp,x,k)
    % EVALMON_DER evaluates the value of partial derivative 
    % of problem in the standard monomial basis.
    %
    %   [value] = EVALMON_DER(coef,supp,x,k) evaluates the problem given by its
    %   internal representation (coefficients and support) in the standard
    %   monomial basis.
    %
    %   Input argument:
    %       coef [double(k,l,m)]: different coefficient/coefficient matrices.
    %       supp [cell(k,n)]: support of the problem.
    %       x [double(1,n)]: point to evaluate.
    %       k : the variable for the partial derivative
    %
    %   Output argument:
    %       value [double]: evaluation.
    %
    %   See also BASISMON.
    
    % MacaulayLab (2023) - Christof Vermeersch.
    % Bor Plestenjak (2026)

    x = x(:).'; % must be a row

    [~, p, q] = size(coef);
    fac = supp(:,k);
    supp(:,k) = max(supp(:,k)-1,0); 
    supp = fac.*prod(x.^supp,2).*ones(1,p,q);
    value = squeeze(sum(coef.*supp,1));
end