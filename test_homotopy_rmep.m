% set to zero to use the same example and experiment with 
% the homotopy part only
constructAB = 1;

if constructAB
    n = 4;
    k = 5;
    rng(1);
    
    A = cell(1,k+1);
    for j = 1:k+1
        A{j} = randn(n+k-1,n);
    end
    
    tic
    [B,Lambda0,X0] = initial_rmep(n,k);
    t2 = toc;
end

lambda02 = remove_duplicates(Lambda0);
nNu = size(Lambda0,1);
disp(['Construction of the initial rmep with ',num2str(nNu),' solutions, ',num2str(length(lambda02)),' distinct : ', num2str(t2), ' s'])
tn = cell(1,nNu);
yn = cell(1,nNu);

lambdaT = [];
resT0 = [];
for j = 1:nNu
    W = B{1};
    for s = 1:k
        W = W + Lambda0(j,s)*B{s+1};
    end
    resT0(j,1) = norm(W*X0(:,j));
end
maxres0 = max(resT0);
disp(['Maximal norm of the initial residual: ',num2str(maxres0)])

opts = [];
opts.display = 1;

if nNu > 4000
    % different parameters when there are many solutions
    opts.innertol = 1e-11;
    opts.maxstepsize = 1e-2;
end

tic
[lambdaT, XT, tn, yn, stat] = homotopy_rmep(A,B,Lambda0,X0,opts);
tTrace = toc;
steps = cellfun(@numel, tn);

fprintf('Homotopy finished with an average of %d steps (max: %d, min %d) in %f s \n',round(mean(steps)),max(steps),min(steps),tTrace)
lambdaT = [];
resT = [];
for j = 1:length(yn)
    lambdaT(j,:) = yn{j}(n+1:end,end).';
    W = lambdaT(j,1)*A{1};
    for s = 1:k
        W = W + lambdaT(j,s+1)*A{s+1};
    end
    resT(j,1) = norm(W*XT(:,j));
end

maxres = max(resT);
lambdaN = lambdaT(:,2:end)./lambdaT(:,1);
lambdaT2 = remove_duplicates(lambdaN);
disp(['Homotopy solver required ', num2str(tTrace+t2), 's to find ',num2str(length(lambdaT2)),'  distinct eigenvalues. Maximal norm of the residual: ',num2str(maxres)])


