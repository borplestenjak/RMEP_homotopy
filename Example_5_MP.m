% Quadratic 2-parameter RMEP solved in multiple precision using Advanpix MCT

rng(1)
A0 = rand(3,2,'mp');
A1 = rand(3,2,'mp');
A2 = rand(3,2,'mp');
A3 = rand(3,2,'mp');
A4 = rand(3,2,'mp');
A5 = rand(3,2,'mp');

A = {A0,A1,A2,A3,A4,A5};
suppA = [0 0;1 0;0 1;2 0;1 1;0 2];
k = 2;
n = 2;

[B,Lambda0,X0] = initial_poly_rmep(2,2,[2 2]);
suppB = [zeros(1,k); 2*eye(k)];
for j=1:numel(B)
    B{j} = numeric_t(B{j},'mp');
end
Lambda0 = numeric_t(Lambda0,'mp');
X0 = numeric_t(X0,'mp');


DA = cell(1,6);
for j=1:6 
    DA{j} = double(A{j});
end

DB = cell(1,3);
for j=1:3 
    DB{j} = double(B{j});
end

opts = [];
opts.display = 1;
fprintf("\nHomotopy in double precision \n----------------------------\n")
[lambdaT, XT] = homotopy_poly_rmep(DA,suppA,DB,suppB,double(Lambda0),double(X0),opts);
lambda = lambdaT(:,2:end)./lambdaT(:,1);

% We compute residuals, geometric multipicites and condition numbers for
% the eigenvalues 
res_gm_cond = []; 
for j = 1:length(lambda)
    res = min(svd(eval_rmep(DA,suppA,lambda(j,:))));
    [gm,s] = condeig_rmep(DA,suppA,lambda(j,:));
    res_gm_cond(j,:) = [res gm s];
end
res_gm_cond_lambda = [res_gm_cond lambda]

opts.innertol = 1e4*eps('mp');
opts.finaltol = 10*eps('mp');

fprintf("\nHomotopy in higher precision \n----------------------------\n")
[lambdaMT, MXT] = homotopy_poly_rmep(A,suppA,B,suppB,Lambda0,X0,opts);
lambdaM = lambdaMT(:,2:end)./lambdaMT(:,1);

% We compute residuals, geometric multipicites and condition numbers for
% the eigenvalues 
res_gm_condM = []; 
for j = 1:length(lambdaM)
    res = min(svd(eval_rmep(A,suppA,lambdaM(j,:))));
    [gm,s] = condeig_rmep(A,suppA,lambdaM(j,:));
    res_gm_condM(j,:) = [res gm s];
end
res_gm_cond_lambdaM = double([res_gm_condM lambdaM])

