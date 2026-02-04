% Linear RMEP with a double eigenvalue that is geometrically simple
%
% This 3x2 linear 2-parameter RMEP is such that normal rank of 
% A0 + lambda1*A1 + lambda2*A2 is 2. It has one geometrically simple
% double eigenvalue (0,0) and a simple eigenvalue (0,1)
%
% Homotopy with default settings repeats the computation of two paths that 
% lead to the double eigenvalue. This could be prevented by setting 
% repeat_opt='simple as then we repeat only paths leading to simple 
% eigenvalues obtained two of more times
 
A0 = [0 0; -1 0; 0 0];
A1 = [1 0;  0 1; 0 0];
A2 = [0 0;  1 0; 0 1];

% Multipareig
fprintf('\nMultipar eig\n------------\n')
lambda = rect_multipareig({A0,A1,A2})

% MacaulayLab
fprintf('\nMacaulayLab \n-----------\n')
lambda1 = rect_multipareig_macaulay({A0, A1, A2})

% Homotopy
fprintf('\nHomotopy 1\n----------\n')
suppA = [0 0; 1 0; 0 1];
A = {A0, A1, A2};
opts = [];
opts.display = 1;
lambda2 = rect_multipareig_homotopy({A0,A1,A2},opts)

% Homotopy, where we set to repeat only simple eigenvalues
fprintf("\nHomotopy 2, repeat_opt = 'simple' \n---------------------------------\n")
opts.repeat_opt = 'simple';
lambda3 = rect_multipareig_homotopy({A0,A1,A2},opts)

% We compute resdiduals, geometric multiplicities and condition numbers of the three
% eigenvalues returned by the homotopy method. 
res_gm_cond3 = []; 
for j = 1:3
    res = min(svd(eval_rmep(A,suppA,lambda3(j,:))));
    [gm,s] = condeig_rmep(A,[],lambda3(j,:));
    res_gm_cond3(j,:) = [res gm s];
end
res_gm_cond_lambda3 = [res_gm_cond3 lambda3]
