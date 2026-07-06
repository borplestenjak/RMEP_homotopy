% Linear RMEP with a finite positive-dimensional solution
%
% This 3x2 linear 2-parameter RMEP is such that the normal rank of 
% A0 + lambda1*A1 + lambda2*A2 is 2. 
% For all points (lambda1,lambda2) such that lambda1+lambda2=1 rank drops
% to 1, so we have a positive-dimensional solution. In addition there is
% one isolated eigenvalue (-4,4), which is a true eigenvalue

% Bor Plestenjak 2026

A0 = [1 3; 
      2 2; 
      3 1];

A1 = [1 2; 
      2 3; 
      3 1];

A2 = [1 1;
      2 2;
      3 0];

% Multipareig finds the isolated point if we flag the problem as singular
fprintf('\nMultiParEig\n-----------\n')
opts = [];
opts.singular = 1;
lambda = rect_multipareig({A0,A1,A2},opts)

% MacaulayLab finds many eigenvalues on the positive-dimensional solution set
fprintf('\nMacaulayLab \n-----------\n')
A = {A0, A1, A2};
suppA = [0 0; 1 0; 0 1];
try
    lambda1 = rect_multipareig_macaulay(A,30,true)
catch ME
    fprintf('Error in macaulaylab: %s \n',ME.message)   
end

% Homotopy finds three eigenvalues, the isolated one and two random ones,
fprintf('\nHomotopy A \n----------\n')
B0A = [2+1i 0; 3-1i 1+1i; 0 2-1i];
B1 =  [1 0; 0 1; 0 0];
B2 =  [0 0; 1 0; 0 1];
BA = {B0A,B1,B2};

[Lambda0,X0] = rect_multipareig(BA)
[lambdaT, XT, ~, ~, stat] = homotopy_rmep(A,BA,Lambda0,X0,opts);
lambdaN1 = lambdaT(:,2:end)./lambdaT(:,1)

fprintf('\nHomotopy B \n----------\n')
B0B = [3+1i 0; 4+1i 5-1i; 0 6-2i];
B1 =  [1 0; 0 1; 0 0];
B2 =  [0 0; 1 0; 0 1];
BB = {B0B,B1,B2};

[Lambda0,X0] = rect_multipareig(BB)
[lambdaT, XT, ~, ~, stat] = homotopy_rmep(A,BB,Lambda0,X0,opts);
lambdaN2 = lambdaT(:,2:end)./lambdaT(:,1)

% We compute resdiduals, geometric multiplicities and condition numbers of the three
% eigenvalues returned by the homotopy method. Al eigenvalues are 
% geometrically simple, but two are highly ill-conditioned, which is a sign 
% that these are not true eigenvalues
res_gm_cond1 = []; 
for j = 1:3
    res = min(svd(eval_rmep(A,suppA,lambdaN1(j,:))));
    [gm,s] = condeig_rmep(A,[],lambdaN1(j,:));
    res_gm_cond1(j,:) = [res gm s];
end
res_gm_cond_lambda1 = [res_gm_cond1 lambdaN1]

res_gm_cond2 = []; 
for j = 1:3
    res = min(svd(eval_rmep(A,suppA,lambdaN2(j,:))));
    [gm,s] = condeig_rmep(A,[],lambdaN2(j,:));
    res_gm_cond2(j,:) = [res gm s];
end
res_gm_cond_lambda2 = [res_gm_cond2 lambdaN2]