% Computation of LTI2 model using cubic two-parameter RMEP by Walsh

addpath examples

N = 14;
rng(1);
y = randn(N,1);

[A,suppA] = lsrwalsh(y,2);

fprintf('\nHomotopy using the Walsh representation\n--------\n')
opts.display = 1;
opts.maxruns = 2;
opts.repeat_opt = 'real';
tic
[points,val,err,cand,simple] = walsh_lti2_homotopy(y); 
t1 = toc
nsolutions = [size(points,1) size(cand,1)]
points
[tmp,pos] = min(val);
solMAC = points(pos,:)
value = tmp

if N<=12
    fprintf('lti2 from MultiParEig\n---------------------\n')
    tic
    [points2,val2,err2,cand2] = lti2(y);
    t2 = toc
    nsolutions2 = [size(points2,1) size(cand2,1)]
    points2
    [tmp,pos] = min(val2);
    solution2 = points2(pos,:)
    value2 = tmp
end

fprintf('MacaulayLab\n---------------------\n')
Amep = mepstruct(A,suppA);
[solMAC,details] = macaulaylab(Amep,posdim=true);
t3 = sum(details.time)
lambdaM = solMAC.num;
indM = find(vecnorm(imag(lambdaM),Inf,2)<1e-8);
nsolutions3 = [size(lambdaM,1) length(indM)]
solution3 = lambdaM(indM,:)


