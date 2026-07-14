function [Res1, Res2, Res3] = compare_rmep_poly_solvers(deg,n,k,mHomo,mMEP,mMac)

% Solve a random polynomial RMEP using three methods and report results 

% Bor Plestenjak 2026

Res1 = [];
Res2 = [];
Res3 = [];

% Options for the homotopy_poly_rmep
opts = [];
opts.display = 1;  
opts.maxruns = 5;
opts.maxangle = 2.5e-1;
opts.maxstepsize = 2.5e-1;
opts.maxinnersteps = 6; %
opts.stepsize = 1e-3;

% Do not change the part below
% ---------------------------------------------------------------------

if deg>2 || k>4
    % MultiParEig solver works only for quadratic polynomials up to
    % degree 4
    mMEP = 0;
end

meja = max([mHomo mMEP mMac]);
m = n + k - 1;
Z0 = zeros(m,n);

for ind = 1:meja

    fprintf('Run %d/%d for deg=%d, n=%d, k=%d\n',ind,meja,deg,n,k)
    % Random full polynomial RMEP
    suppA = monomials(deg,k);
    suppA = suppA(:,k:-1:1);
    nA = size(suppA,1);
    
    rng(ind); 
    A = cell(1,nA);
    for j = 1:nA
        A{j} = randn(n+k-1,n);
    end

    if ind<=mHomo
        % use homotopy to solve the problem

        % start parallel pool if not active
        p = gcp('nocreate');
        if isempty(p)
            parpool
        end

        rng(100+ind)
        tic
        [B,Lambda0,X0] = initial_poly_rmep(n,k,deg*ones(1,k));
        suppB = [zeros(1,k); deg*eye(k)];
        nNu = size(Lambda0,1);
        init_t = toc;

        fprintf('Construction of the initial rmep with %d solutions in %7.3e sec\n',nNu, init_t);
        tn = cell(1,nNu);
        yn = cell(1,nNu);

        % Homotopy method for PRMEP
        tic
        [lambdaT, XT, tn, yn, stat] = homotopy_poly_rmep(A,suppA,B,suppB,Lambda0,X0,opts);
        tTrace = toc;
        steps = cellfun(@numel, tn);
        lambda = lambdaT(:,2:end)./lambdaT(:,1); % conversion to affine eigenvalues
        % fprintf('Homotopy finished with an average of %d steps (max: %d, min %d) in %f s \n',round(mean(steps)),max(steps),min(steps),tTrace)
        
        resT = 0;
        for j = 1:size(lambda,1)
            W = eval_rmep(A,suppA,lambda(j,:),Z0);
            if n>1
                resT(j,1) = norm(W*XT(:,j))/norm(W);
            else
                resT(j,1) = norm(W*XT(:,j));
            end
        end
        
        Res1.full_time(ind) = tTrace + init_t;
        maxres = max(resT);
        Res1.maxres(ind) = max(resT);
        Res1.meanres(ind) = median(resT);
        lambda2 = remove_duplicates(lambda);
        Res1.distinct(ind) = length(lambda2);
        Res1.init_time(ind) = init_t;
        Res1.avgsteps(ind) = mean(steps);
        if stat.run>1
            Res1.repeat_paths(ind) = sum(stat.paths(2:end));
        else
            Res1.repeat_paths(ind) = 0;
        end
        disp(['Homotopy    solver required ', num2str(tTrace+init_t), 's to find ',num2str(length(lambda2)),'  distinct eigenvalues. Maximal norm of the residual: ',num2str(maxres)])
    end

    if ind<=mMEP
        tic; 
        try
            [lambda,X] = rect_quad_multipareig(A); 
        catch ME
             fprintf('Error in multipareig: %s \n',ME.message)   
             lambda = [];
             X = [];
        end
            
        tM = toc;

        resT = 0;
        for j = 1:size(lambda,1)
            W = eval_rmep(A,suppA,lambda(j,:),Z0);
            if n>1
                resT(j,1) = norm(W*X(:,j))/norm(W);
            else
                resT(j,1) = norm(W*X(:,j));
            end
        end
        maxres = max(resT);
        Res2.full_time(ind) = tM;
        Res2.maxres(ind) = max(resT);
        Res2.meanres(ind) = median(resT);
        lambda2 = remove_duplicates(lambda);
        Res2.distinct(ind) = length(lambda2);
        disp(['Multipareig solver required ', num2str(tM), 's to find ',num2str(length(lambda2)),'  distinct eigenvalues. Maximal norm of the residual: ',num2str(maxres)])
    end

    if ind<=mMac
        tic; 
        [lambda,X,details] = poly_rect_multipareig_macaulay(A,suppA); 
        tM = toc;

        resT = 0;
        for j = 1:size(lambda,1)
            W = eval_rmep(A,suppA,lambda(j,:),Z0);
            if n>1
                resT(j,1) = norm(W*X(:,j))/norm(W);
            else
                resT(j,1) = norm(W*X(:,j));
            end
        end
        maxres = max(resT);
        Res3.full_time(ind) = tM;
        Res3.maxres(ind) = max(resT);
        Res3.meanres(ind) = median(resT);
        lambda2 = remove_duplicates(lambda);
        Res3.distinct(ind) = length(lambda2);
        disp(['MacaulayLab solver required ', num2str(tM), 's to find ',num2str(length(lambda2)),'  distinct eigenvalues. Maximal norm of the residual: ',num2str(maxres)])
    end
    
end

if mHomo>0
    Res1.repeat_avg = mean(Res1.repeat_paths);    
    Res1.avgsteps_avg = mean(Res1.avgsteps);    
    Res1.init_time_avg = mean(Res1.init_time);
    Res1.full_time_avg = mean(Res1.full_time);
    Res1.maxres_avg = mean(Res1.maxres);    
    Res1.meanres_avg = mean(Res1.meanres);
end

if mMEP>0
    Res2.full_time_avg = mean(Res2.full_time);
    Res2.maxres_avg = mean(Res2.maxres);    
    Res2.meanres_avg = mean(Res2.meanres);
end

if mMac>0
    Res3.full_time_avg = mean(Res3.full_time);
    Res3.maxres_avg = mean(Res3.maxres);    
    Res3.meanres_avg = mean(Res3.meanres);
end
