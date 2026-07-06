% Linear RMEP with an infinite eigenvalue
%
% This 3x2 linear 2-parameter RMEP A0 + lambda1*A1 + lambda2*A2
% is such that the rank of A2 is 1, therefore the homogeneous version
% has eigenvalue (0,0,1), which corresponds to an infinite eigenvalue. 

% Bor Plestenjak 2026

rng(1)
A0 = rand(3,2);
A1 = rand(3,2);
A2 = 1*ones(3,2);

A0 = [ 1 -1
      -2  1
      -1  2];

A1 = [ 2 -1
      -2  1
       4 -1];

B0A = [2 0;
       3 1;
       0 2];
B1 = [1 0;0 1;0 0];
B2 = [0 0;1 0;0 1];

[lambda0,X0] = rect_multipareig({B0A,B1,B2})

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

% Homotopy finds both finite eigenvalues and another one at infinity
opts = [];
opts.display = 4;
fprintf('\nHomotopy \n----------\n')
[lambda2,X2,lambdaT2,XT2] = rect_multipareig_homotopy(A,opts);
lambda2
lambdaT2

