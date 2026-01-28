% Finds solution of LTI(2) using homotoy for polynomial RMEPs

sampleN = 8; % size of the sample 
constructAB = 1; % set to zero to use the same example and experiment with the homotopy part only

if constructAB

    rng(1)
    y = randn(sampleN,1);

    [M00,M10,M01,M20,M11,M02] = LTI2_mat(y);
    [n1,n] = size(M00);

    k = 2;  % number of parameters
    deg = 2;  % degree
    
    suppA = [0 0; 1 0; 0 1; 2 0; 1 1; 0 2];
    nA = size(suppA,1);  

    rng(1);
    A = {M00,M10,M01,M20,M11,M02};

    tic
    [B,Lambda0,X0] = initial_poly_rmep(n,k,[2 2]);
    suppB = [0 0;2 0;0 2];
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

% parameters have to be adjusted from case to case
opts = [];
opts.display = 1;
opts.maxruns = 2; 
opts.maxstepsize = 1e-1;
opts.maxangle = 1e-1; 
opts.stepsize = 1e-12;
opts.repeat_opt = 'real'; % heuristic method where we do not repeat complex multiple eigenvalzues that have mathing number of conjugate values

% for larger samples we usually get correct number of real eigenvalues and
% very close number of finite eigenvalues. This might be improved if we
% reduce maxangle, maxstepsize, and stepsize
if sampleN>9
    % tighter options increases precison but also computational time
    opts.maxstepsize = 2.5e-2;
    opts.maxangle = 5e-2; % reducing this increases precison and computational time
end

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

maxT = cellfun(@max, tn);
indT = find(maxT>1-1e-12);

maxres = max(resT);
lambdaN = lambdaT(:,2:end)./lambdaT(:,1);
lambdaT2 = remove_duplicates(lambdaN);
disp(['Homotopy solver required ', num2str(tTrace+t2), 's to find ',num2str(length(lambdaT2)),'  distinct eigenvalues. Maximal norm of the residual: ',num2str(maxres)])

% extract finite solutions
indF = find(abs(lambdaT(indT,1))>1e-2); % this setting might not be good for all finite solutions
ind = indT(indF);
cand = lambdaN(ind,:);
maxresfin = max(resT(ind));
disp(['Found ', num2str(length(cand)),'  finite eigenvalues. Maximal norm of the residual: ',num2str(maxresfin)])

% extract real solutions
ind2 = find(abs(imag(cand(:,1)))+abs(imag(cand(:,2)))<1e-6);
disp(['Found ', num2str(length(ind2)),'  real eigenvalues.  Maximal norm of the residual: ',num2str(max(resT(ind(ind2))))])
real_solutions = real(cand(ind2,:))
