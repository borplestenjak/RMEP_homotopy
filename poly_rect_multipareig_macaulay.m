function [lambda,X,details] = poly_rect_multipareig_macaulay(A,suppA,dend,posdim)

if nargin<33
    dend = max(size(A{1},1)+5,100);
end

if nargin<4
    posdim = 0;
end

k = numel(A) - 1;
Amep = mepstruct(A,suppA);

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
            W = eval_rmep(A,suppA,lambda(i,:));
            X(:,i) = min_sing_vec(W,0); % we have to use SVD as matrices are rectangular
        end
    end
catch ME
    fprintf('Error in macaulaylab: %s \n',ME.message)   
end

