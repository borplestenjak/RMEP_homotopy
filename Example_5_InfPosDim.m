% Quadratic 2-parameter RMEP with positive-dimensional solution in infinity
% 
% We take 3x2 quadratic 2-parameter RMEP of the form 
% such that for all the homogoeneous version is singular for all nonzero
% (0,a,b) such that a^2+ab+b^2=0.
%
% The problem has 10 finite eigenvalues (a generic QRMEP of this size has 
% 12 solutions) 
% 
% If we take initial QRMEP of form B0+lambda1^2*B1+lambda2^2*B2 then
% 2 paths go to infinite solutions (we get different points, 
% depending on the initial RMEP B).
%
% By computing geometric multiplicities and condition number of the 
% computed 12 eigenvalues, we clearly see which ones are true eigenvalues
% and which are from paths leading to infinity

rng(1)
A0 = rand(3,2);
A1 = rand(3,2);
A2 = rand(3,2);
A3 = rand(3,2);
A4 = [A3(:,1) rand(3,1)];
A5 = [A3(:,1) rand(3,1)];

A = {A0,A1,A2,A3,A4,A5};
suppA = [0 0; 1 0; 0 1; 2 0; 1 1;0 2];

% Multipareig finds 6 true eigenvalues
fprintf('\nMultiParEig\n-----------\n')
lambda = rect_quad_multipareig(A)

% In two runs with the initial QRMEP of form B0+lambda1^2*B1+lambda2^2*B2
% we get 12 solutions that include 6 true eigenvalues
opts = [];
opts.display = 2;

% Homotopy
fprintf('\nHomotopy 1 using initial RMEP B0 + lambda1^2*B1 + lambda2^2*B2 \n--------------------------------------------------------------\n')
rng(2)
[B1,Lambda10,X10] = initial_poly_rmep(2,2,[2 2]);
suppB1 = [0 0;2 0;0 2];
[lambdaT1, XT1] = homotopy_poly_rmep(A,suppA,B1,suppB1,Lambda10,X10,opts);
lambda1 = lambdaT1(:,[2 3])./lambdaT1(:,1)

fprintf('\nHomotopy 2 using initial RMEP B0 + lambda1^2*B1 + lambda2^2*B2 \n--------------------------------------------------------------\n')
rng(3)
[B2,Lambda20,X20] = initial_poly_rmep(2,2,[2 2]);
suppB2 = [0 0;2 0;0 2];
[lambdaT2, XT2] = homotopy_poly_rmep(A,suppA,B2,suppB2,Lambda20,X20,opts);
lambda2 = lambdaT2(:,[2 3])./lambdaT2(:,1)

% We compute resdiduals, geometric multiplicities and condition numbers of the three
% eigenvalues returned by the homotopy method. 
res_gm_cond2 = []; 
for j = 1:12
    res = min(svd(eval_rmep(A,suppA,lambda2(j,:))));
    [gm,s] = condeig_rmep(A,suppA,lambda2(j,:));
    res_gm_cond2(j,:) = [res gm s];
end
res_gm_cond_lambda2 = [res_gm_cond2 lambda2]

% MacaulayLab 
fprintf("\nMacaulayLab using posdim=true \n---------------------------------\n")
Amep = mepstruct(A,suppA);

options = [];
options.posdim = true;
try
    solution = macaulaylab(Amep,posdim=true);
    lambdaM  = solution.num
    resM = [];
    for j = 1:size(lambdaM,1)
        W = eval_rmep(A,suppA,lambdaM(j,:));
        resM(j,1) = min(svd(W));
    end
    maxresfinM = max(resM);
    disp(['Found ', num2str(length(lambdaM)),'  finite eigenvalues. Maximal norm of the residual: ',num2str(maxresfinM)])
catch ME
    fprintf('Error in macaulaylab: %s \n',ME.message)   
end
