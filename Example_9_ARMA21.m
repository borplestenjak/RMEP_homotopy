% Computation of ARMA(2,1) and comparison to solution obtained by MultiParEig

N = 5;
rng(1);
y = randn(N,1);

fprintf('\nHomotopy\n--------\n')
tic
[points,val,err,cand,simple] = arma21_homotopy(y); 
t1 = toc
nsolutions = [size(points,1) size(cand,1)]
points
[tmp,pos] = min(val);
solution = points(pos,:)
value = tmp


fprintf('lti2 from MultiParEig\n---------------------\n')
tic
[points2,val2,err2,cand2] = arma21(y);
t2 = toc
nsolutions2 = [size(points2,1) size(cand2,1)]
points2
[tmp,pos] = min(val2);
solution2 = points2(pos,:)
value2 = tmp
 
% find differences between finite solutions found
[Da,Db] = compare_sets(cand,cand2,1e-6); Da, Db
