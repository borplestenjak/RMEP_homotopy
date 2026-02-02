% set to zero to use the same example and experiment with 
% the homotopy part only
constructAB = 1;

if constructAB
    N = 4;
    na = 1; 
    nc = 1;

    [A, suppA] = arma(randn(N,1), na, nc);
    n = size(A{1}, 2);
    k = na + nc;
    deg = max(sum(suppA, 2));

    tic
    [B,Lambda0,X0] = initial_poly_rmep(n,k,deg*ones(1,k));
    suppB = [zeros(1,k); deg*eye(k)];
    t2 = toc;
    m = n + k - 1;
end

lambda02 = remove_duplicates(Lambda0);
nNu = size(Lambda0,1);
disp(['Construction of the initial rmep with ',num2str(nNu),' solutions, ',num2str(length(lambda02)),' distinct : ', num2str(t2), ' s'])
tn = cell(1,nNu);
yn = cell(1,nNu);
Z0 = zeros(m,n);

lambdaT = [];
resT0 = [];
for j = 1:nNu
    W = eval_rmep(B,suppB,lambda02(j,:),Z0);
    resT0(j,1) = norm(W*X0(:,j));
end
maxres0 = max(resT0);
disp(['Maximal norm of the initial residual: ',num2str(maxres0)])

opts = [];
opts.maxruns = 1;
opts.display = 1;

if nNu > 4000
     % different parameters when there are many solutions
     opts.innertol = 1e-11;
     opts.maxstepsize = 1e-2;
end

tic
[lambdaT, XT, tn, yn, stat] = homotopy_poly_rmep(A,suppA,B,suppB,Lambda0,X0,opts);
tTrace = toc;
steps = cellfun(@numel, tn);

fprintf('Homotopy finished with an average of %d steps (max: %d, min %d) in %f s \n',round(mean(steps)),max(steps),min(steps),tTrace)

% homogeneous representations of A
homosuppA = [deg-sum(suppA,2) suppA];

for j = 1:length(yn)
    W = eval_rmep(A,homosuppA,lambdaT(j,:),Z0);
    resT(j,1) = norm(W*XT(:,j));
end

maxres = max(resT);
lambdaN = lambdaT(:,2:end)./lambdaT(:,1);
lambdaT2 = remove_duplicates(lambdaN);
disp(['Homotopy solver required ', num2str(tTrace+t2), 's to find ',num2str(length(lambdaT2)),'  distinct eigenvalues. Maximal norm of the residual: ',num2str(maxres)])
