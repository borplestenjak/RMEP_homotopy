function [lambda,X,details] = rect_multipareig_macaulay(A,dend,posdim)

if nargin<2
    dend = max(size(A{1},1)+5,30);
end

if nargin<3
    posdim = 0;
end

k = numel(A) - 1;
supp = [zeros(1,k); eye(k)];
Amep = mepstruct(A,supp);

lambda = [];
X = [];
details = [];

try
    if posdim
        [solution,details] = macaulaylab(Amep,dend,posdim=true); 
    else
        [solution,details] = macaulaylab(Amep,dend); 
    end
    lambda = solution.num;
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
catch ME
    fprintf('Error in macaulaylab: %s \n',ME.message)   
end

