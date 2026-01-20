function [lambda,X,tcell,ycell,stat] = homotopy_poly_rmep(coefA,suppA_r,coefB,suppB_r,Lambda0,X0,opts)

% HOMOTOPY_POLY_RMEP   Solve a polynomial rectangular multiparameter 
% eigenvalue problem using homotopy
%
% [lambda,X,tn,yn,stat] = homotopy_rmep(coefA,suppA,coefB,suppB,Lambda0,X0,opts) returns 
% eigenvalues and eigenvectors of a polynomial rectangular multiparameter eigenvalue 
% problem (RMEP) 
%
% (A{1} + lambda(1) A{2} x + ... + lambda(k) A{k+1}) x = 0 
%
% using initial RMEP
%
% (B{1} + lambda(1) B{2} x + ... + lambda(k) B{k+1}) x = 0 
%
% with known eigenvalues Lambda0 and eigenvectors X0
%
% Input:
%   - coefA,suppA: polynomial RMEP, all matrices have to be 
%         rectangular matrices of the same size (n+k-1) x n
%   - coefB,suppB: initial polynomial RMEP, all matrices have to be 
%         rectangular matrices of the same size (n+k-1) x n
%   - Lambda0: rows are eigenvalues of the RMEP B (matrix of size m x k),
%         where m = nchoosek(n+k-1,k)
%   - X0: columns are eigenvectors of the RMEP B (matrix of size n x m)
%   - opts : options 
%
% Options in opts:
%   - stepsize (1e-4): initial stepsize for the homotopy
%   - minstepsize (1e-10): minimal stepsize for the homotopy
%   - maxstepsize (2.5e-2): maximal stepsize for the homotopy
%   - maxsteps (20000): maximum number of steps in the homotopy
%   - maxinnersteps (6): maximum number of steps for Newton's corrector
%   - innertol (3e-13): relative tolerance for inner homotopy steps
%   - finaltol (1e-15): relative tolerance for final homotopy step
%   - maxruns (4): maximal number of attempts to follow curves
%   - display (0): display all iterations (2), statistics of each run (1), nothing (0) 
%   - dist_tol (1e-6): relative tolerance for checking if the solutions are distinct
%   - abort_stepsize (1e-20): tolerance to abort path following if the stepsize becomes to small
%
% Output:
%   - lambda : matrix m x k, each row is an eigenvalue
%   - X : matrix n x m with right eigenvectors
%   - tcell: cell with all time steps 
%   - ycell: cell with all intermediate steps, first n elements represent eigenvector,
%         last k elements represent eigenvalue
%   - stat: several statistics

% Bor Plestenjak, 2025, 2026

if nargin < 7, opts = []; end

if isfield(opts,'stepsize'),      stepsize = opts.stepsize;           else,  stepsize = 1e-4;      end
if isfield(opts,'minstepsize'),   minstepsize = opts.minstepsize;     else,  minstepsize = 1e-14;  end
if isfield(opts,'maxstepsize'),   maxstepsize = opts.maxstepsize;     else,  maxstepsize = 2.5e-2;   end
if isfield(opts,'maxsteps'),      maxsteps = opts.maxsteps;           else,  maxsteps = 20000;      end
if isfield(opts,'maxinnersteps'), maxinnersteps = opts.maxinnersteps; else,  maxinnersteps = 6;    end
if isfield(opts,'innertol'),      innertol = opts.innertol;           else,  innertol = sqrt(10)*1e-13;     end
if isfield(opts,'finaltol'),      finaltol = opts.finaltol;           else,  finaltol = 1e-15;     end
if isfield(opts,'maxruns'),       maxruns = opts.maxruns;             else,  maxruns = 6;          end
if isfield(opts,'display'),       display = opts.display;             else,  display = 0;          end
if isfield(opts,'dist_tol'),      dist_tol = opts.dist_tol;           else,  dist_tol = 1e-6;      end
if isfield(opts,'abort_stepsize'), abort_stepsize = opts.abort_stepsize;     else,  abort_stepsize = 1e-20;  end

kpar = size(suppA_r,2);
[m,n] = size(squeeze(coefA(1,:,:)));
neig = size(Lambda0,1); % number of curves to follow
tcell = cell(1,neig);
ycell = cell(1,neig);
indices = 1:neig; % indices of curves that we have to follow

stat = [];
stat.fail = 0;
stat.follow_ind = indices;
stat.indices = [];
stat.run = 0;

converged = zeros(1,neig);
minh = ones(1,neig);
run = 0;

deg = max(sum(suppA_r,2)); % degree of the RMEP
% homogeneous representations of RMEPS A and B
suppA = [deg-sum(suppA_r,2) suppA_r];
suppB = [deg-sum(suppB_r,2) suppB_r];

while run<maxruns && sum(converged)<neig
    run = run + 1; 
    follow_ind = find(converged==0);
    minh(follow_ind) = 1;

    % we follow paths in parallel
    for ind = 1:neig
        if converged(ind)==0
            Ma = zeros(m+2,n+kpar+1);
            b = zeros(m+2,1);
            Ma0 = zeros(m+2,n+kpar+1);
            b0 = zeros(m+2,1);
            C = cell(1,kpar+1);   
            % we follow homogeneous RMEP
            lambda0 = [1; Lambda0(ind,:).'];
            lambda0 = lambda0/norm(lambda0);
            t = 0;
            y = [X0(:,ind); lambda0]; % initial vector for the homotopy
            tn = t;
            yn = y;
            step = 0;
            h = stepsize;
            abort_path = 0;
            while (t<1-1e-12) && (step<maxsteps) && (abort_path==0)
                step = step + 1;
                tol = innertol;
                if t + h >= 1
                    h = 1 - t;
                    tol = finaltol;
                end

                W = (1-t)*evalmon(coefB,suppB,y(n+1:end)) + t*evalmon(coefA,suppA,y(n+1:end));
                Wt = evalmon(coefA,suppA,y(n+1:end)) - evalmon(coefB,suppB,y(n+1:end));
            
                Ma0(1:m,1:n) = W; 
                for j = 1:kpar+1
                    Ma0(1:m,n+j) = ((1-t)*evalmon_der(coefB,suppB,y(n+1:end),j) + t*evalmon_der(coefA,suppA,y(n+1:end),j))*y(1:n);
                end
                Ma0(m+1,1:n) = 2*y(1:n)';
                Ma0(m+2,n+1:end) = 2*y(n+1:end)';
                b0(1:m) = -Wt*y(1:n);
            
                warning off
                tang = Ma0\b0; % tangent vector (with t'=1)
                warning on

                noconv = 1;
                corrector_tol = tol;
                while noconv && (abort_path==0)
            
                    % compute predictor
                    yp = y + h*tang;
                    tp = t + h;
                    % compute corrector by Newton's method
                    
                    % initial residual
                    W = (1-tp)*evalmon(coefB,suppB,yp(n+1:end)) + tp*evalmon(coefA,suppA,yp(n+1:end));
                    rhs = [W*yp(1:n); yp(1:n)'*yp(1:n)-1; yp(n+1:end)'*yp(n+1:end)-1];
                    innerstep = 0;
                    
                    % Newton correction
                    while (norm(rhs)>corrector_tol*norm(W)) && (innerstep<maxinnersteps)
                        innerstep = innerstep + 1;
                        Ma(1:m,1:n) = W;
                        for j = 1:kpar+1
                            Ma(1:m,n+j) = ((1-tp)*evalmon_der(coefB,suppB,yp(n+1:end),j) + tp*evalmon_der(coefA,suppA,yp(n+1:end),j))*yp(1:n);
                        end
                        Ma(m+1,1:n) = 2*yp(1:n)';
                        Ma(m+2,n+1:end) = 2*yp(n+1:end)';
                        warning off
                        dx = -Ma\rhs;
                        warning on
                        yp = yp + dx; 
                        W = (1-tp)*evalmon(coefB,suppB,yp(n+1:end)) + tp*evalmon(coefA,suppA,yp(n+1:end));
                        rhs = [W*yp(1:n); yp(1:n)'*yp(1:n)-1; yp(n+1:end)'*yp(n+1:end)-1];
                    end
                    if norm(rhs)<corrector_tol*norm(W)
                        t = tp;
                        y = yp;
                        tn = [tn t];
                        yn = [yn y];
                        noconv = 0;
                        % if convergence is fast, we increase h 
                        if innerstep < maxinnersteps-1
                            h = min(2*h,maxstepsize);
                        end
                    else
                        % no convergence, we reduce h and repeat from the
                        % last good point on the path
                        h = h/2;
                        minh(ind) = min(minh(ind),h);
                        if h<abort_stepsize
                            abort_path = 1;
                        end
                    end
                end
            end
            tcell{ind} = tn;
            ycell{ind} = yn;
            X(:,ind) = yn(1:n,end);
            lambda(ind,:) = yn(n+1:end,end).';
            if display>1
                disp([ind step tn(end) lambda(ind,:)])
            end
            if t>1-1e-12
                converged(ind) = 1;
            end
        end
    end    
    lambdaN = lambda(:,2:end)./lambda(:,1);
    indices = find_duplicates(lambdaN,dist_tol); % if we find duplicates, we repeat the computations for these indices using smaller stepsize
    fail = 0;
    for j = 1:neig
        if (tcell{j}(end)<1-1e-12) && (length(find(indices==j))==0)
            fail = fail + 1;
            if length(find(indices==j))==0
                indices = sort([indices j]);
            end
        end
    end
    stat.fail = fail;
    stat.follow_ind = follow_ind;
    stat.indices = indices;
    stat.run = run;
    if display
        nconv = sum(converged(follow_ind));
        steps = cellfun(@numel, tcell);
        steps = steps(follow_ind);
        minstepused = min(minh(follow_ind));
        fprintf('Run %d, converged %d out of %d paths, fail %d, to repeat %d, average %d steps (max: %d, min %d) minstep %6.2e, initial stepsize %6.2e, maxstepsize %6.2e, tol %6.2e, maxNewton %d\n',...
            run, nconv, length(follow_ind), fail, length(indices), round(mean(steps)),max(steps),min(steps),minstepused,stepsize,maxstepsize,innertol,maxinnersteps)
    end
    if ~isempty(indices)
        converged(indices)=0;
        % tighter criteria for the repeated curves
        stepsize = max(stepsize/100,minstepsize);
        maxstepsize = maxstepsize/3;
        innertol = max(innertol/sqrt(10),sqrt(10)*1e-15);
        maxinnersteps = max(maxinnersteps-1,4);
    end
end


