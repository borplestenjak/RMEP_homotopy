% Finds solution of ARMA(1,1) using homotoy for polynomial RMEPs

sampleN = 6; % size of the sample
constructAB = 1; % set to zero to use the same example and experiment with the homotopy part only

if constructAB

    rng(2)
    y = randn(sampleN,1);

    [M00,M10,M01,M02] = ARMA11_mat(y);
    [n1,n] = size(M00);

    k = 2;  % number of parameters
    deg = 2;  % degree
    
    suppA = [0 0; 1 0; 0 1; 0 2];
    nA = size(suppA,1);

    rng(1);
    A = cell(1,nA);
    A{1} = M00;
    A{2} = M10;
    A{3} = M01;
    A{4} = M02;

    tic
    [B,Lambda0,X0] = initial_poly_rmep(n,k,[1 2]);
    suppB = [0 0;1 0;0 2];
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
opts.display = 1;
opts.maxruns = 1; % for ARMA RMEP almost all points qualify for repeated computations, so we do 1 run only
opts.use_pinv = 1;
opts.maxstepsize = 1e-1;
opts.maxangle = 1e-1;
opts.stepsize = 1e-12;
opts.repeat_opt = 'real'; % heuristic method where we do not repeat complex multiple eigenvalzues that have mathing number of conjugate values

tic
[lambdaT, XT, tn, yn, stat] = homotopy_poly_rmep(A,suppA,B,suppB,Lambda0,X0,opts);
tTrace = toc;
steps = cellfun(@numel, tn);

fprintf('Homotopy finished with an average of %d steps (max: %d, min %d) in %f s \n',round(mean(steps)),max(steps),min(steps),tTrace)

% homogeneous representations of A
homosuppA = [deg-sum(suppA,2) suppA];

resT = [];
for j = 1:length(yn)
    W = eval_rmep(A,homosuppA,lambdaT(j,:),Z0);
    resT(j,1) = norm(W*XT(:,j));
end

maxres = max(resT);
lambdaN = lambdaT(:,2:end)./lambdaT(:,1);
lambdaT2 = remove_duplicates(lambdaN);
disp(['Homotopy solver required ', num2str(tTrace+t2), 's to find ',num2str(length(lambdaT2)),'  distinct eigenvalues. Maximal norm of the residual: ',num2str(maxres)])

% extract finite solutions
ind = find(abs(lambdaT(:,1))>0.1); 
cand = lambdaN(ind,:);
maxresfin = max(resT(ind));
disp(['Found ', num2str(length(cand)),'  finite eigenvalues. Maximal norm of the residual: ',num2str(maxresfin)])

% extract real solutions
ind2 = find(abs(imag(cand(:,1)))+abs(imag(cand(:,2)))<1e-6);
disp(['Found ', num2str(length(ind2)),'  real eigenvalues.  Maximal norm of the residual: ',num2str(max(resT(ind(ind2))))])
real_solutions = real(cand(ind2,:))
