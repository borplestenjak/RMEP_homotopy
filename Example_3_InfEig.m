% Quadratic 2-parameter RMEP with only one nonzero monomial of degree 2
% 
% We take 3x2 quadratic 2-parameter RMEP of the form 
% A0 + lambda1*A1 + lambda2*A2 + lambda1^2*A3
%
% Since none of the qudratic monomials containing lambda2 is present, there
% are only 6 eigenvalues (a generic QRMEP of this size has 12 solutions). A
% homogeneous version has infinite eigenvalue (0,0,1) of multiplicity 6.
% 
% If we take initial QRMEP of form B0+lambda1^2*B1+lambda2^2*B2 then
% 6 paths lead to 6 eigenvalues and the remaining 6 paths go to infinity 
% (we get different points, depending on the initial RMEP B)
%
% If we take initial QRMEP of form B0+lambda1^2*B1+lambda2*B2, then
% 6 paths lead to 6 eigenvalues
%
% By computing geometric multiplicities and condition number of the 
% computed 12 eigenvalues, we clearly see which ones are true eigenvalues
% and which are from paths leading to infinity

rng(1)
A0 = rand(3,2);
A1 = rand(3,2);
A2 = rand(3,2);
A3 = rand(3,2);
A4 = zeros(3,2);
A5 = zeros(3,2);

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

% If we use the initial QRMEP of form B0+lambda1^2*B1+lambda2*B2
% we get only 6 true eigenvalues
fprintf('\nHomotopy 3 using initial RMEP B0 + lambda1^2*B1 + lambda2*B2 \n------------------------------------------------------------\n')
rng(4)
[B3,Lambda30,X30] = initial_poly_rmep(2,2,[2 1]);
suppB3 = [0 0;2 0;0 1];
[lambdaT3, XT3] = homotopy_poly_rmep(A,suppA,B3,suppB3,Lambda30,X30,opts);
lambda3 = lambdaT3(:,[2 3])./lambdaT3(:,1)

% We compute residuals, geometric multipicites and condition numbers for 12
% eigenvalues obtained in the second run with the initial QRMEP of 
% form B0+lambda1^2*B1+lambda2^2*B2
res_gm_cond2 = []; 
for j = 1:length(lambda2)
    res = min(svd(eval_rmep(A,suppA,lambda2(j,:))));
    [gm,s] = condeig_rmep(A,suppA,lambda2(j,:));
    res_gm_cond2(j,:) = [res gm s];
end
res_gm_cond_lambda2 = [res_gm_cond2 lambda2]
