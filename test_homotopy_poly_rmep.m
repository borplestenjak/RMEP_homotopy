% set to zero to use the same example and experiment with 
% the homotopy part only
constructAB = 1;

if constructAB
    n = 10;
    k = 2;
    deg = [2 2];
    rng(1);
    coefA = [];
    coefB = [];
    suppA = [];
    suppB = [];
    
    A = cell(1,6);
    for j = 1:6
        A{j} = randn(n+k-1,n);
        coefA(j,:,:) = A{j}; 
    end
    suppA = [0 0;1 0;0 1;2 0;1 1;0 2]; % random 2-parameter QRMEP
    
    tic
    [B,Lambda0,X0] = initial_poly_rmep(n,k,deg);
    for j = 1:k+1;
        coefB(j,:,:) = B{j};
    end
    suppB = [0 0;2 0; 0 2];
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
    W = evalmon(coefB,suppB,lambda02(j,:));
    resT0(j,1) = norm(W*X0(:,j));
end
maxres0 = max(resT0);
disp(['Maximal norm of the initial residual: ',num2str(maxres0)])

opts = [];
opts.display = 1;

% if nNu > 4000
%     % different parameters when there are many solutions
%     opts.innertol = 1e-11;
%     opts.maxstepsize = 1e-2;
% end

tic
[lambdaT, XT, tn, yn, stat] = homotopy_poly_rmep(coefA,suppA,coefB,suppB,Lambda0,X0,opts);
tTrace = toc;
steps = cellfun(@numel, tn);

fprintf('Homotopy finished with an average of %d steps (max: %d, min %d) in %f s \n',round(mean(steps)),max(steps),min(steps),tTrace)

maxdeg = max(sum(suppA,2)); % degree of the RMEP
% homogeneous representations of RMEPS A and B
homosuppA = [maxdeg-sum(suppA,2) suppA];

for j = 1:length(yn)
    W = evalmon(coefA,homosuppA,lambdaT(j,:));
    resT(j,1) = norm(W*XT(:,j));
end

maxres = max(resT);
lambdaN = lambdaT(:,2:end)./lambdaT(:,1);
lambdaT2 = remove_duplicates(lambdaN);
disp(['Homotopy solver required ', num2str(tTrace+t2), 's to find ',num2str(length(lambdaT2)),'  distinct eigenvalues. Maximal norm of the residual: ',num2str(maxres)])
 

