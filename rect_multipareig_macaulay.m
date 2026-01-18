function [lambda,X] = rect_multipareig_macaulay(A)

k = size(A,2) - 1;
Amep = mepstruct(A,1,k);
lambda = macaulaylab(Amep); 
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
