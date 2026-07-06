% Computation of ARMA(1,1) and comparison to solution obtained by MultiParEig

N = 8;
rng(1);
y = randn(N,1);

fprintf('\nHomotopy\n--------\n')
tic
[points,val,err,cand,simple] = arma11_homotopy(y); 
t1 = toc
nsolutions = [size(points,1) size(cand,1)]
format long e
points
[tmp,pos] = min(val);
solution1 = points(pos,:)
value1 = tmp
format short e


fprintf('arma11 from MultiParEig\n---------------------\n')
opta = [];
opta.showrank = 1;
tic
[points2,val2,err2,cand2] = arma11(y,opta);
t2 = toc
nsolutions2 = [size(points2,1) size(cand2,1)]
format long e
points2
[tmp,pos] = min(val2);
solution2 = points2(pos,:)
value2 = tmp
format short e

[M00,M10,M01,M02] = ARMA11_mat(y);
mep = {M00,M10,M01,M02};
suppA = [0 0;1 0;0 1;0 2];
delta = sqrt(eps);

if N<=7
    fprintf('macaulaylab from MacaulayLab\n---------------------\n')
    Amep = mepstruct(mep, suppA);
    
    try
        tic
        solution = macaulaylab(Amep,posdim=true);
        t3 = toc
        cand3 = solution.num;
        ind3 = find(abs(imag(cand3(:,1)))+abs(imag(cand3(:,2)))<delta)
        points3 = cand3(ind3,:)
        nsolutions3 = size(cand3,1)
    catch ME
        fprintf('Error in macaulaylab: %s \n',ME.message)   
    end
end

% find differences between finite solutions found
[Da,Db] = compare_sets(cand,cand2,1e-5);
Da_length = size(Da,1), 
Db_length = size(Db,1)

% Plot points
for param = 1:2
    figure
    plot(real(points2(:,param)),imag(points2(:,param)),'om','MarkerSize',20)
    hold on
    plot(real(cand2(:,param)),imag(cand2(:,param)),'.g','MarkerSize',40)
    plot(real(cand2(:,param)),imag(cand2(:,param)),'.g','MarkerSize',40)
    plot(real(cand(:,param)),imag(cand(:,param)),'.k','MarkerSize',25)
    if N<=7
        plot(real(cand3(:,param)),imag(cand3(:,param)),'.r','MarkerSize',12)
    end
end

res1 = residuals_rmep(mep,suppA,cand);
res1sol = residuals_rmep(mep,suppA,solution1);
res2 = residuals_rmep(mep,suppA,cand2);
res2sol = residuals_rmep(mep,suppA,solution2);
if N<=7
    res3 = residuals_rmep(mep,suppA,cand3);
    res3sol = residuals_rmep(mep,suppA,points3);
else
    res3 = Inf;
    res3sol = Inf;
end

max_residual = [max(res1) max(res2) max(res3)]
avg_residual = [mean(res1) mean(res2) mean(res3)]

res1sol
res2sol
res3sol
