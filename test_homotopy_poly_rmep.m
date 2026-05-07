% Solve a random polynomial RMEP using homotopy

% Change the folowing parameters 
% ---------------------------------------------------------------------
n = 6;    % number of columns in matrices A
k = 4;    % number of parameters
deg = 3;  % degree of the polynomial RMEP
rg1 = 1;  % random generator for the construction of A
rg2 = 2;  % random generator for the construction of the initial problem B

% Options for the homotopy_poly_rmep
opts = [];
opts.display = 1;
opts.maxruns = 3;

% Do not change the part below
% ---------------------------------------------------------------------

% Random full polynomial RMEP
suppA = monomials(deg,k);
nA = size(suppA,1);

rng(rg1);
A = cell(1,nA);
for j = 1:nA
    A{j} = randn(n+k-1,n);
end
    
% Construction of the initial problem
tic
rng(rg2)
[B,Lambda0,X0] = initial_poly_rmep(n,k,deg*ones(1,k));
suppB = [zeros(1,k); deg*eye(k)];
t2 = toc;
m = n + k - 1;

nNu = size(Lambda0,1);
disp(['Construction of the initial rmep with ',num2str(nNu),' solutions: ', num2str(t2), ' s'])
tn = cell(1,nNu);
yn = cell(1,nNu);
Z0 = zeros(m,n);

lambdaT = [];
resT0 = [];
for j = 1:nNu 
    W = eval_rmep(B,suppB,Lambda0(j,:),Z0);
    resT0(j,1) = norm(W*X0(:,j));
end
maxres0 = max(resT0);
disp(['Maximal norm of the initial residual: ',num2str(maxres0)])

% Homotopy method for PRMEP
tic
[lambdaT, XT, tn, yn, stat] = homotopy_poly_rmep(A,suppA,B,suppB,Lambda0,X0,opts);
tTrace = toc;
steps = cellfun(@numel, tn);
lambda = lambdaT(:,2:end)./lambdaT(:,1); % conversion to affine eigenvalues
fprintf('Homotopy finished with an average of %d steps (max: %d, min %d) in %f s \n',round(mean(steps)),max(steps),min(steps),tTrace)

resT = [];
for j = 1:length(yn)
    W = eval_rmep(A,suppA,lambda(j,:),Z0);
    resT(j,1) = norm(W*XT(:,j));
end

maxres = max(resT);
lambda2 = remove_duplicates(lambda);
disp(['Homotopy solver required ', num2str(tTrace+t2), 's to find ',num2str(length(lambda2)),'  distinct eigenvalues. Maximal norm of the residual: ',num2str(maxres)])
