% Computation of H2-norm reduced-order model
% It runs the identification problem for different data lengths.

addpath examples

% Set parameters:
maxm = 3;
nmin = 7;
nmax = 9;
delta = 1e-5;

timehomotopy = zeros(maxm-1,nmax);
timemultipareig = zeros(maxm-1,nmax);
timemacaulaylab = zeros(maxm-1,nmax);
errhomotopy = zeros(maxm-1,nmax);
errmultipareig = zeros(maxm-1,nmax);
errmacaulaylab = zeros(maxm-1,nmax);

for n = nmin:nmax
 fprintf('n = %d\n',n)
 disp('=================================')

% Run the different experiments:
for m = 2:maxm
    fprintf("Run for reduced order %d \n", m)
        rng(1)
        a = randn(n,1);
        b = randn(n,1);
        [mat, supp] = h2sisored3(a,b,m);
        mep = mepstruct(mat, supp);
        
        fprintf("Homotopy \n")
        opts.display = 1;
        opts.maxruns = 2;
        opts.repeat_opt = 'real';
        opts.abort_inf = 1e-4;  
        opts.maxstepsize = 2e-1;
        opts.maxangle = 1e-1;   
        opts.stepsize = 1e-8;   
        tstart = tic;
        [points1, ~, lambdaT,~] = poly_rect_multipareig_homotopy(mat,supp,opts);
        t1 = toc(tstart)
        timehomotopy(m-1,n) = t1;
        ind = find(vecnorm(imag(points1),Inf,2) < delta);
        points1 = real(points1(ind, :));
        nsol1 = length(points1)
        errhomotopy(m-1, n) = residuals(points1, mep);
        
        if m<=4
            fprintf("MultiParEig (standard approach) \n")
            options.singular = 1;
            options.showrank = 1;
            tstart = tic;
            points2 = rect_quad_multipareig(mat, options);
            points2 = points2(:,m:-1:1); % supp ordering is different in MultiParEig
            t2 = toc(tstart)
            timemultipareig(m-1,n) = t2;
            ind = find(vecnorm(imag(points2),Inf,2) < delta);
            points2 = real(points2(ind, :));
            nsol2 = length(points2)
            errmultipareig(m-1, n) = residuals(points2, mep);
        end

        if n<25
            fprintf("MacaulayLab (standard approach) \n")
            tstart = tic;
            points3 = macaulaylab(mep, posdim = true);
            t3 = toc(tstart)
            timemacaulaylab(m-1, n) = t3;
            points3 = num(points3);
            ind = find(vecnorm(imag(points3),Inf,2) < delta);
            points3 = real(points3(ind, :));
            nsol3 = length(points3)
            errmacaulaylab(m-1, n) = residuals(points3, mep);
        end
    end
end

timeh = timehomotopy(:,nmin:nmax)
timeMEP = timemultipareig(:,nmin:nmax)
timeMac = timemacaulaylab(:,nmin:nmax)
errh = errhomotopy(:,nmin:nmax)
errMEP = errmultipareig(:,nmin:nmax)
errMac = errhomotopy(:,nmin:nmax)