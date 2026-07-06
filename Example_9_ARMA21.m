% Computation of ARMA(2,1) and comparison to solution obtained by MultiParEig

N = 6;
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


fprintf('arma21 from MultiParEig\n---------------------\n')
tic
opts_arma = [];
opts_arma.showrank = 1;
[points2,val2,err2,cand2] = arma21(y,opts_arma);
t2 = toc
nsolutions2 = [size(points2,1) size(cand2,1)]
points2
[tmp,pos] = min(val2);
solution2 = points2(pos,:)
value2 = tmp
 
% find differences between finite solutions found
[Da,Db] = compare_sets(cand,cand2,1e-6); 
Da_length = size(Da,1), 
Db_length = size(Db,1)
