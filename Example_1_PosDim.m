% Linear RMEP with a finite positive-dimensional solution
%
% This 3x2 linear 2-parameter RMEP is such that the normal rank of 
% A0 + lambda1*A1 + lambda2*A2 is 2. 
% For all points (lambda1,lambda2) such that lambda1+lambda2=1 rank drops
% to 1, so we have a positive-dimensional solution. In addition there is
% one isolated eigenvalue (0.9580,-0.6894), which is a true eigenvalue

% Bor Plestenjak 2026

rng(1)
A0 = rand(3,2);
A1 = [-A0(:,1) rand(3,1)];
A2 = [-A0(:,1) rand(3,1)];

% Multipareig finds the isolated point if we flag the problem as singular
fprintf('\nMultiParEig\n-----------\n')
opts = [];
opts.singular = 1;
lambda = rect_multipareig({A0,A1,A2},opts)

% MacaulayLab does not find the gap, even if we use the flag posdim=true
fprintf('\nMacaulayLab \n-----------\n')
A = {A0, A1, A2};
suppA = [0 0; 1 0; 0 1];
try
    lambda1 = rect_multipareig_macaulay(A,30,true)
catch ME
    fprintf('Error in macaulaylab: %s \n',ME.message)   
end

% Homotopy finds three eigenvalues, the isolated one and two random ones,
% if we run the homotopy again with a different initial problem, the
% solutions have only the true eigenvalue in common. 
opts = [];
opts.display = 1;
opts.maxruns = 1;
fprintf('\nHomotopy 1 \n----------\n')
lambda2 = rect_multipareig_homotopy(A,opts)
fprintf('\nHomotopy 2 \n----------\n')
lambda3 = rect_multipareig_homotopy(A,opts)

% We compute resdiduals, geometric multiplicities and condition numbers of the three
% eigenvalues returned by the homotopy method. Al eigenvalues are 
% geometrically simple, but two are highly ill-conditioned, which is a sign 
% that these are not true eigenvalues
res_gm_cond3 = []; 
for j = 1:3
    res = min(svd(eval_rmep(A,suppA,lambda3(j,:))));
    [gm,s] = condeig_rmep(A,[],lambda3(j,:));
    res_gm_cond3(j,:) = [res gm s];
end
res_gm_cond_lambda3 = [res_gm_cond3 lambda3]