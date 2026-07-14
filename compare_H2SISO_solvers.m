function [Res1, Res2, Res3] = compare_H2SISO_solvers(m,n,mHomo,mMEP,mMac)

% Solve a H2SISO problem and compare results

Res1 = [];
Res2 = [];
Res3 = [];

% Options for the homotopy_poly_rmep
opts = [];
opts.display = 1;
opts.maxruns = 3;
opts.abort_inf = 1e-4;  
opts.maxstepsize = 2e-1;
opts.maxangle = 1e-1;   
opts.stepsize = 1e-8;   

delta = 1e-5;
% Do not change the part below
% ---------------------------------------------------------------------

meja = max([mHomo mMEP mMac]);
% m = n + k - 1;

for run = 1:meja

    fprintf('Run %d/%d for m=%d, n=%d\n',run,meja,m,n)
    % construction of a random linear RMEP
    rng(run)
    a = randn(n,1);
    b = randn(n,1);
    [mat, supp] = h2sisored3(a,b,m);
    mep = mepstruct(mat, supp);
 
    if run<=mHomo
        % use homotopy to solve the problem

        % start parallel pool if not active
        p = gcp('nocreate');
        if isempty(p)
            parpool
        end

        % Homotopy method for PRMEP
        tstart = tic;
        points1 = poly_rect_multipareig_homotopy(mat,supp,opts);
        tTrace = toc(tstart)
        indh = find(vecnorm(imag(points1),Inf,2) < delta);
        points1 = real(points1(indh, :));
        nsol1 = length(points1)
        Res1.full_time(run) = tTrace;
        Res1.maxres(run) = residuals(points1, mep);
        disp(['Homotopy    solver required ', num2str(tTrace), 's'])
    end

    if run<=mMEP
        options = [];
        options.singular = 1;
        options.showrank = 1;
        options.rankeps = 1e-9;
        tstart = tic;
        points2 = rect_quad_multipareig(mat, options);
        t2 = toc(tstart)
        Res2.full_time(run) = t2;
        if length(points2)>0
            points2 = points2(:,m:-1:1); % supp ordering is different in MultiParEig
            indh = find(vecnorm(imag(points2),Inf,2) < delta);
            points2 = real(points2(indh, :));
            nsol2 = length(points2)
            Res2.maxres(run) = residuals(points2, mep);
        else
            Res2.maxres(run) = 0;
        end
        disp(['Multipareig solver required ', num2str(t2), 's'])
    end
        
    if run<=mMac
        tstart = tic;
        points3 = macaulaylab(mep, posdim = true);
        t3 = toc(tstart)
        points3 = num(points3);
        indh = find(vecnorm(imag(points3),Inf,2) < delta);
        points3 = real(points3(indh, :));
        nsol3 = length(points3)
        Res3.full_time(run) = t3;
        Res3.maxres(run) = residuals(points3, mep);
        disp(['MacaulayLab solver required ', num2str(t3), 's'])
    end
    
end

if mHomo>0
    Res1.full_time_avg = mean(Res1.full_time);
    Res1.maxres_avg = mean(Res1.maxres);    
end

if mMEP>0
    Res2.full_time_avg = mean(Res2.full_time);
    Res2.maxres_avg = mean(Res2.maxres);    
end

if mMac>0
    Res3.full_time_avg = mean(Res3.full_time);
    Res3.maxres_avg = mean(Res3.maxres);    
end
