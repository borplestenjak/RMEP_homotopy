function [value] = evalmon(coef,supp,x)
    % EVALMON evaluates the problem in the standard monomial basis.
    %
    %   [value] = EVALMON(coef,supp,x) evaluates the problem given by its
    %   internal representation (coefficients and support) in the standard
    %   monomial basis.
    %
    %   Input argument:
    %       coef [double(k,l,m)]: different coefficient/coefficient matrices.
    %       supp [cell(k,n)]: support of the problem.
    %       x [double(1,n)]: point to evaluate.
    %
    %   Output argument:
    %       value [double]: evaluation.
    %
    %   See also BASISMON.
    
    % MacaulayLab (2023) - Christof Vermeersch.
    
    x = x(:).'; % must be a row
    
    [~, p, q] = size(coef);
    supp = prod(x.^supp,2).*ones(1,p,q);
    value = squeeze(sum(coef.*supp,1));
end