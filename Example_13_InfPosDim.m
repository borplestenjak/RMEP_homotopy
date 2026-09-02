% Linear 3-parameter RMEP 
% 
% We take linear 3-parameter RMEP of the form
% A0 + lambda1 * A1 + lambda2 * A2 + lambda3 * A3

A0=[2 2 3
    7 2 3
    2 2 9
    2 4 8
    1 4 6];

A1=[2 6 2
    3 4 3
    0 0 4
    0 0 1
    0 0 7];

A2=[5 4 7
   10 2 3
    0 0 6
    0 0 9
    0 0 6];

A3=[4 6 1
    3 8 3
    0 0 4
    0 0 7
    0 0 10];

A = {A0,A1,A2,A3};
k = 3;
n = 2;

% Multipareig finds the isolated point if we flag the problem as singular
fprintf('\nMultiParEig\n-----------\n')
opts = [];
opts.singular = 1;
lambda = rect_multipareig(A,opts)

% MacaulayLab 
fprintf('\nMacaulayLab \n-----------\n')
suppA = [0 0 0; 1 0 0; 0 1 0;0 0 1];
try
    lambda1 = rect_multipareig_macaulay(A,50,true)
catch ME
    fprintf('Error in macaulaylab: %s \n',ME.message)   
end

% Homotopy finds three eigenvalues, the isolated one and two random ones,
% if we run the homotopy again with a different initial problem, the
% solutions have only the true eigenvalue in common. 
opts = [];
opts.display = 3;
opts.maxruns = 1;
fprintf('\nHomotopy 1 \n----------\n')
[lambda2,X2,lambdaT2,XT2] = rect_multipareig_homotopy(A,opts);
fprintf('\nHomotopy 2 \n----------\n')
[lambda3,X3,lambdaT3,XT3] = rect_multipareig_homotopy(A,opts);

% We compute residuals, geometric multiplicities and condition numbers of the three
% eigenvalues returned by the homotopy method. Al eigenvalues are 
% geometrically simple, but two are highly ill-conditioned, which is a sign 
% that these are not true eigenvalues
res_gm_cond2 = []; 
for j = 1:size(lambdaT2,1)
    [gm,s] = condeig_rmep(A,[],lambdaT2(j,:).',false);
    res_gm_cond2(j,:) = [gm s];
end
res_gm_cond_lambda2 = [res_gm_cond2 lambdaT2]