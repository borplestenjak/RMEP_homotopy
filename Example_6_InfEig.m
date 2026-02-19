% Linear RMEP with an infinite eigenvalue
%
% This 3x2 linear 2-parameter RMEP A0 + lambda1*A1 + lambda2*A2
% is such that the rank of A2 is 1, therefore the homogeneous version
% has eigenvalue (0,0,1), which corresponds to an infinite eigenvalue. 

% Bor Plestenjak 2026

rng(1)
A0 = rand(3,2);
A1 = rand(3,2);
A2 = ones(3,2);

% Multipareig finds the finite eigenvalues if we flag the problem as singular
fprintf('\nMultiParEig\n-----------\n')
opts = [];
opts.singular = 1;
lambda = rect_multipareig({A0,A1,A2},opts)

% MacaulayLab
fprintf('\nMacaulayLab \n-----------\n')
A = {A0, A1, A2};
suppA = [0 0; 1 0; 0 1];
options = [];
options.posdim = true;
try
    lambda1 = rect_multipareig_macaulay(A,30,true)
catch ME
    fprintf('Error in macaulaylab: %s \n',ME.message)   
end

% Homotopy finds both finite eigenvalues, the isolated one and two random ones,
% if we run the homotopy again with a different initial problem, the
% solutions have only the true eigenvalue in common
opts = [];
opts.display = 1;
fprintf('\nHomotopy \n----------\n')
lambda2 = rect_multipareig_homotopy(A,opts)

% We compute resdiduals, geometric multiplicities and condition numbers of the three
% eigenvalues returned by the homotopy method. Al eigenvalues are 
% geometrically simple, but two are highly ill-conditioned, which is a sign 
% that these are not true eigenvalues
res_gm_cond2 = []; 
for j = 1:size(lambda2,1)
    res = min(svd(eval_rmep(A,suppA,lambda2(j,:))));
    [gm,s] = condeig_rmep(A,[],lambda2(j,:));
    res_gm_cond2(j,:) = [res gm s];
end
res_gm_cond_lambda2 = [res_gm_cond2 lambda2]