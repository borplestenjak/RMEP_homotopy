function [lambda,X] = rect_multipareig_macaulay(A,dend,opts)

if nargin<2
    dend = 30;
end

if nargin<3
    opts = [];
end

k = numel(A) - 1;
Amep = mepstruct(A,1,k);
if isempty(opts)
    lambda = macaulaylab(Amep,dend); 
else
    lambda = macaulaylab(Amep,dend,opts); 
end

n = size(A{1},2);

if nargout>1
    X = zeros(n,size(lambda,1));
    for i = 1:size(lambda,1)
        TMP = A{1};
        for j = 1:k
            TMP = TMP + lambda(i,j)*A{j+1};
        end
        X(:,i) = min_sing_vec(TMP,0); % we have to use SVD as matrices are rectangular
    end
end
