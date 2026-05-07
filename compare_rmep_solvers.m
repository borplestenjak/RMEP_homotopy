function [Res1, Res2, Res3] = compare_rmep_solvers(n,k,mHomo,mMEP,mMac)

% Solve a random linear RMEP using three methods and report results 

% Bor Plestenjak 2026

Res1 = [];
Res2 = [];
Res3 = [];

% Options for the homotopy_poly_rmep
opts = [];
opts.display = 1;  
opts.maxruns = 4;
opts.maxangle = 2.5e-1;
opts.maxstepsize = 2.5e-1;
opts.maxinnersteps = 6; %

% Do not change the part below
% ---------------------------------------------------------------------

meja = max([mHomo mMEP mMac]);
m = n + k - 1;
Z0 = zeros(m,n);

for ind = 1:meja

    fprintf('Run %d/%d for n=%d, k=%d\n',ind,meja,n,k)
    % construction of a random linear RMEP
    rng(ind); 
    A = cell(1,k+1);
    for j = 1:k+1
        A{j} = randn(n+k-1,n);
    end   
 
    if ind<=mHomo
        % use homotopy to solve the problem

        % start parallel pool if not active
        p = gcp('nocreate');
        if isempty(p)
            parpool
        end

        max_init_cond = Inf;
        init_run = 0;
        rng(100+ind)
        test_cond = 0;

        tic
        if test_cond
            while init_run<5 && max_init_cond>1e28
    
                init_run = init_run +1 ; 
                [B,Lambda0,X0] = initial_poly_rmep(n,k,ones(1,k));
                nNu = size(Lambda0,1);
                suppB = [];
                % test initial RMEP
                init_cond = [];
                for j = 1:nNu 
                    [gm,s] = condeig_rmep(B,[],Lambda0(j,:));
                    init_cond(j,1) = s;
                end
                max_init_cond = max(init_cond)
            end
        else
            [B,Lambda0,X0] = initial_poly_rmep(n,k,ones(1,k));
            nNu = size(Lambda0,1);
        end
        init_t = toc;

        if test_cond
           fprintf('Construction of the initial rmep with %d solutions in %d runs, %7.3e sec, max_cond %7.3e\n',nNu, init_run, init_t,max_init_cond);
        else
            fprintf('Construction of the initial rmep with %d solutions in %7.3e sec\n',nNu, init_t);
        end
        % disp(['Construction of the initial rmep with ',num2str(nNu),' solutions: ', num2str(init_t), ' s, max_cond:', num2str(max_init_cond) ])
        tn = cell(1,nNu);
        yn = cell(1,nNu);

        % % test initial RMEP
        % lambdaT = [];
        % resT0 = [];
        % for j = 1:nNu 
        %     W = eval_rmep(B,[],Lambda0(j,:),Z0);
        %     resT0(j,1) = norm(W*X0(:,j));
        % end
        % maxres0 = max(resT0);
        % disp(['Maximal norm of the initial residual: ',num2str(maxres0), ', maximal cond: ',num2str(max_init_cond)])

        % Homotopy method for PRMEP
        tic
        [lambdaT, XT, tn, yn, stat] = homotopy_rmep(A,B,Lambda0,X0,opts);
        tTrace = toc;
        steps = cellfun(@numel, tn);
        lambda = lambdaT(:,2:end)./lambdaT(:,1); % conversion to affine eigenvalues
        % fprintf('Homotopy finished with an average of %d steps (max: %d, min %d) in %f s \n',round(mean(steps)),max(steps),min(steps),tTrace)
        
        resT = 0;
        for j = 1:size(lambda,1)
            W = eval_rmep(A,[],lambda(j,:),Z0);
            resT(j,1) = norm(W*XT(:,j))/norm(W);
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
            [lambda,X] = rect_multipareig(A); 
        catch ME
             fprintf('Error in multipareig: %s \n',ME.message)   
             lambda = [];
             X = [];
        end
            
        tM = toc;

        resT = 0;
        for j = 1:size(lambda,1)
            W = eval_rmep(A,[],lambda(j,:),Z0);
            resT(j,1) = norm(W*X(:,j))/norm(W);
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
        [lambda,X,details] = rect_multipareig_macaulay(A); 
        tM = toc;

        resT = 0;
        for j = 1:size(lambda,1)
            W = eval_rmep(A,[],lambda(j,:),Z0);
            resT(j,1) = norm(W*X(:,j))/norm(W);
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
