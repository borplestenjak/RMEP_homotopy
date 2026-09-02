% Cubic 3-parameter RMEP with several zero monomials 
% 
% We take cubic 3-parameter RMEP of the form
% A0 + lambda1*A1 + lambda1^2*A2 + lambda2^2*A3 + lambda3^3*A4
%
% Due to the structure the RMEP has only 48 eigenvalues (a generic cubic
% 3-parameter RMEP of this size has 108 solutions). The homogeneous version
% has positive-dimensional solution (0,lambda1,lambda2,0) for all nonzero
% (lambda1,lambda2). 
%
% If we start homotopy from B0 + lambda1^3*B1 + lambda2^3*B2 + lambda3^3*B3,
% we follow 108 paths to get 48 eigenvalues
% If we start homotopy from B0 + lambda1^2*B1 + lambda2^2*B2 + lambda3^3*B3,
% we follow 48 paths to get 48 eigenvalues

% Bor Plestenjak 2026

rng(1)
m = 2;
A0 = randn(m+2,m);
A1 = randn(m+2,m);
A2 = randn(m+2,m);
A3 = randn(m+2,m);
A4 = randn(m+2,m);

A = {A0,A1,A2,A3,A4};
suppA = [0 0 0;1 0 0;2 0 0;0 2 0;0 0 3];
k = 3;
n = 2;

fprintf('\nHomotopy 1 using initial RMEP B0 + lambda1^3*B1 + lambda2^3*B2 + lambda3^3*B3\n')
fprintf('-----------------------------------------------------------------------------\n')
[B,Lambda0,X0] = initial_poly_rmep(m,3,[3 3 3]);
suppB = [zeros(1,k); 3*eye(k)];

opts = [];
opts.display = 2;
opts.maxruns = 2;
[lambdaT,XT] = homotopy_poly_rmep(A,suppA,B,suppB,Lambda0,X0,opts);

% extract finite solutions
ind = find(abs(lambdaT(:,1))>1e-4); 
lambda = lambdaT(ind,2:end)./lambdaT(ind,1)
X = XT(:,ind);

Z0 = zeros(m+2,m);
res = [];
for j = 1:length(lambda)
    W = eval_rmep(A,suppA,lambda(j,:),Z0);
    res(j,1) = norm(W*X(:,j));
end

maxresfin = max(res);
disp(['Found ', num2str(length(lambda)),'  finite eigenvalues. Maximal norm of the residual: ',num2str(maxresfin)])

fprintf('\nHomotopy 2 using initial RMEP B0 + lambda1*B1 + lambda2^2*B2 + lambda3^3*B3\n')
fprintf('---------------------------------------------------------------------------\n')
[lambda2,X2] = poly_rect_multipareig_homotopy(A,suppA,opts);

lambda2
res2 = [];
for j = 1:length(lambda2)
    W = eval_rmep(A,suppA,lambda2(j,:),Z0);
    res2(j,1) = norm(W*X2(:,j));
end

maxresfin2 = max(res2);
disp(['Found ', num2str(length(lambda2)),'  finite eigenvalues. Maximal norm of the residual: ',num2str(maxresfin2)])

% MacaulayLab 
fprintf("\nMacaulayLab using posdim=true \n---------------------------------\n")
Amep = mepstruct(A,suppA);

try
    solution = macaulaylab(Amep,30,posdim=true);
    lambdaM = solution.num
    resM = [];
    for j = 1:length(lambda2)
        W = eval_rmep(A,suppA,lambdaM(j,:),Z0);
        resM(j,1) = min(svd(W));
    end
    maxresfinM = max(resM);
    disp(['Found ', num2str(length(lambdaM)),'  finite eigenvalues. Maximal norm of the residual: ',num2str(maxresfinM)])
catch ME
    fprintf('Error in macaulaylab: %s \n',ME.message)   
end



