% Computation of LTI(K) model using cubic K-parameter RMEP by Walsh

addpath examples

K = 3;
N = 10;
rng(1);
y = randn(N,1);

[A,suppA] = lsrwalsh(y,K);

fprintf('\nHomotopy using the Walsh representation\n--------\n')
opts.display = 1;
opts.maxruns = 1;
tic
[cand,X,lambdaT,XT] = poly_rect_multipareig_homotopy(A,suppA,opts);
t1 = toc
ind = find(vecnorm(imag(cand),Inf,2)<1e-8);
nsolutions = [size(cand,1) length(ind)]
points = real(cand(ind,:))

fprintf('MacaulayLab\n---------------------\n')
Amep = mepstruct(A,suppA);
[solMAC,details] = macaulaylab(Amep,posdim=true);
t2 = sum(details.time)
lambdaM = solMAC.num;
indM = find(vecnorm(imag(lambdaM),Inf,2)<1e-8);
nsolutions2 = [size(lambdaM,1) length(indM)]
points2 = real(lambdaM(indM,:))
