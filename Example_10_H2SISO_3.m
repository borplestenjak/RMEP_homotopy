% Computation of H2siso3 model and comparison to solution obtained by MultiParEig

addpath examples

N = 7;
rng(1);
a = randn(N,1);
b = randn(N,1);
[A,suppA] = h2sisored3(a,b,2);

fprintf('\nHomotopy\n--------\n')
opts.display = 1;
opts.maxruns = 2;
opts.repeat_opt = 'real';
tic
[lambda,X,lambdaT,XT] = poly_rect_multipareig_homotopy(A,suppA,opts); toc, size(lambda)
t1 = toc
realind = find(vecnorm(imag(lambda),Inf,2)<1e-8);
nsolutions = [size(lambda,1) length(realind)]
real_solutions = real(lambda(realind,:))

fprintf('lti2 from MultiParEig\n---------------------\n')
tic
[mu,X] = rect_quad_multipareig(A); 
t2 = toc
mu = mu(:,[2 1]); % switch because of different ordering of monomials 
realind2 = find(vecnorm(imag(mu),Inf,2)<1e-8);
nsolutions = [size(mu,1) length(realind2)]
real_solutions2 = real(mu(realind2,:))
 
% find differences between finite solutions found
[Da,Db] = compare_sets(lambda,mu,1e-6); Da, Db
