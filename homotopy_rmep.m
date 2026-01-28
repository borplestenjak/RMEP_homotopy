function [lambda,X,tcell,ycell,stat] = homotopy_rmep(A,B,Lambda0,X0,opts)

% HOMOTOPY_RMEP   Solve a linear rectangular multiparameter eigenvalue 
% problem using homotopy
%
% [lambda,X,tn,yn,stat] = homotopy_rmep(A,B,Lambda0,X0,opts) returns 
% homogeneous eigenvalues and eigenvectors of a rectangular multiparameter 
% eigenvalue problem (RMEP)
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
%   - A,B : cell arrays of size k+1, each cell is an (n+k-1) x n matrix 
%   - Lambda0: rows are eigenvalues of the RMEP B (matrix of size m x (k+1)),
%         where m = nchoosek(n+k-1,k)
%   - X0: columns are eigenvectors of the RMEP B (matrix of size n x m)
%   - opts : options 
%
% Options in opts:
%   - stepsize (1e-8): initial stepsize for the homotopy
%   - minstepsize (2*eps): minimal relative stepsize for the homotopy
%   - maxstepsize (2.5e-2): maximal stepsize for the homotopy
%   - maxsteps (2000): maximum number of steps in the homotopy
%   - maxinnersteps (6): maximum number of steps for Newton's corrector
%   - innertol (3e-13): relative tolerance for inner homotopy steps
%   - finaltol (1e-15): relative tolerance for final homotopy step
%   - maxruns (4): maximal number of attempts to follow curves
%   - display (0): display all iterations (3),  exceptions (2), statistics of each run (1), no display (0) 
%   - dist_tol (1e-6): relative tolerance for checking if the solutions are distinct
%   - fail_stepsize (1e-30): tolerance to abort path following if the stepsize becomes to small
%   - abort_inf (1e-10): tolerance on abs(lambda_0) to identify an infinite eigenvalue
%   - maxangle (1e-1): maximal angle (subspace) between eigenvectors and eigenvalues
%   - repeat_opt ('all'): how to select paths to repeat, 'all': all duplicates, 
%       'real' all real duplicates and complex duplicates without conjugate match
%   - use_pinv (0): solve Jacobian with pinv when the system is close to singular
%
% Output:
%   - lambda : matrix m x (k+1), each row is a homogeneous eigenvalue
%   - X : matrix n x m with right eigenvectors
%   - tcell: cell with all time steps 
%   - ycell: cell with all intermediate steps, first n elements represent eigenvector,
%         last k+1 elements represent eigenvalue
%   - stat: several statistics

% Bor Plestenjak, 2025, 2026

% Validate number of input parameters
narginchk(4, 5);
if nargin < 5, opts = []; end

if isfield(opts,'stepsize'),      stepsize = opts.stepsize;           else,  stepsize = 1e-8;            end
if isfield(opts,'minstepsize'),   minstepsize = opts.minstepsize;     else,  minstepsize = 2*eps;        end
if isfield(opts,'maxstepsize'),   maxstepsize = opts.maxstepsize;     else,  maxstepsize = 2.5e-2;       end
if isfield(opts,'maxsteps'),      maxsteps = opts.maxsteps;           else,  maxsteps = 2000;            end
if isfield(opts,'maxinnersteps'), maxinnersteps = opts.maxinnersteps; else,  maxinnersteps = 6;          end
if isfield(opts,'innertol'),      innertol = opts.innertol;           else,  innertol = 3e-13;           end
if isfield(opts,'finaltol'),      finaltol = opts.finaltol;           else,  finaltol = 1e-15;           end
if isfield(opts,'maxruns'),       maxruns = opts.maxruns;             else,  maxruns = 4;                end
if isfield(opts,'display'),       disp_level = opts.display;          else,  disp_level = 0;             end
if isfield(opts,'dist_tol'),      dist_tol = opts.dist_tol;           else,  dist_tol = 1e-6;            end
if isfield(opts,'fail_stepsize'), fail_stepsize = opts.fail_stepsize; else,  fail_stepsize = 1e-30;      end
if isfield(opts,'abort_inf'),     abort_inf = opts.abort_inf;         else,  abort_inf = 1e-10;          end
if isfield(opts,'maxangle'),      maxangle = opts.maxangle;           else,  maxangle = 1e-1;            end
if isfield(opts,'repeat_opt'),    repeat_opt = opts.repeat_opt;       else,  repeat_opt = 'all';         end
if isfield(opts,'use_pinv'),      use_pinv = opts.use_pinv;           else,  use_pinv = 0;               end

kpar = length(A) - 1;   % number of parameters
[m,n] = size(A{1});
neig = size(Lambda0,1); % number of paths to follow
tcell = cell(1,neig);
ycell = cell(1,neig);
rep_indices = 1:neig; % indices of paths that we have to follow

stat = [];
stat.fail = 0;
stat.follow_ind = rep_indices;
stat.indices = [];
stat.run = 0;

converged = zeros(1,neig);
minh = ones(1,neig);
run = 0;

while run<maxruns && length(find(converged==1))<neig
    run = run + 1; 
    follow_ind = find(converged==0);
    minh(follow_ind) = 1;

    % we follow paths in parallel
    parfor ind = 1:neig
        if converged(ind)==0
            % initialization of data for path following
            Ma = zeros(m+2,n+kpar+1);
            Ma0 = zeros(m+2,n+kpar+1);
            b0 = zeros(m+2,1);
            C = cell(1,kpar+1);   
            % we follow homogeneous RMEP, we convert initial eigenvalue to homogeneous form
            lambda0 = [1; Lambda0(ind,:).'];
            lambda0 = lambda0/norm(lambda0);
            t = 0;
            y = [X0(:,ind); lambda0]; % initial vector for the homotopy
            tn = t;
            yn = y;
            step = 0;
            h = stepsize;
            abort_path = 0;

            % path following
            while (t<1-1e-15) && (step<maxsteps) && (abort_path==0)

                step = step + 1;
                tol = innertol;
                if t + h >= 1
                    h = 1 - t;
                    tol = finaltol; % we increase tolerance for the last step
                end

                W = 0;
                Wt = 0;
                for j = 1:kpar+1
                    C{j} = (1-t)*B{j} + t*A{j}; 
                    W = W + y(n+j)*C{j};
                    Wt = Wt + y(n+j)*(A{j} - B{j}); 
                end
            
                Ma0(1:m,1:n) = W; 
                for j = 1:kpar+1
                    Ma0(1:m,n+j) = C{j}*y(1:n);
                end
                Ma0(m+1,1:n) = 2*y(1:n)';
                Ma0(m+2,n+1:end) = 2*y(n+1:end)';
                b0(1:m) = -Wt*y(1:n);
            
                % tangent vector (with t'=1) for the predictor
                [tang, rcnd] = linsolve(Ma0,b0); 
                if use_pinv && rcnd<1e-14
                    tang = pinv(Ma0)*b0;
                end

                noconv = 1;
                while noconv && (abort_path==0)
            
                    % compute predictor by Euler's method
                    yp = y + h*tang;
                    tp = t + h;

                    % initial residual
                    W = 0;
                    for j = 1:kpar+1
                        C{j} = (1-tp)*B{j} + tp*A{j}; 
                        W = W + yp(n+j)*C{j};
                    end
                    rhs = [W*yp(1:n); yp(1:n)'*yp(1:n)-1; yp(n+1:end)'*yp(n+1:end)-1];
                    innerstep = 0;
                    
                    % compute corrector by Newton's method
                    while (norm(rhs)>tol*norm(W)) && (innerstep<maxinnersteps)
                        innerstep = innerstep + 1;
                        Ma(1:m,1:n) = W;
                        for j = 1:kpar+1
                            Ma(1:m,n+j) =  C{j}*yp(1:n);
                        end
                        Ma(m+1,1:n) = 2*yp(1:n)';
                        Ma(m+2,n+1:end) = 2*yp(n+1:end)';
                        % solve the system with the Jacobian matrix Ma
                        [dx, rcnd] = linsolve(Ma,-rhs);
                        if use_pinv && rcnd<1e-14
                            dx = -(pinv(Ma)*rhs);
                        end
                        yp = yp + dx; 
                        W = 0;
                        for j = 1:kpar+1
                            C{j} = (1-tp)*B{j} + tp*A{j}; 
                            W = W + yp(n+j)*C{j};
                        end
                        rhs = [W*yp(1:n); yp(1:n)'*yp(1:n)-1; yp(n+1:end)'*yp(n+1:end)-1];
                    end

                    angle_vec = norm(y(1:n)-(y(1:n)'*yp(1:n))*yp(1:n));
                    angle_eig = norm(y(n+1:end)-(y(n+1:end)'*yp(n+1:end))*yp(n+1:end));
                    if norm(rhs)<tol*(1+norm(W)) && angle_vec<maxangle && angle_eig<maxangle
                        % Newton method converged, angles are not too large, we take step forward
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
                        % no convergence, we repeat predictor-corrector step with smaller h
                        h = h/2;
                        tol = innertol; 
                        minh(ind) = min(minh(ind),h);
                        % if h is too small, we abort this path
                        if h<(t*minstepsize + fail_stepsize)
                           if disp_level==2 
                                fprintf('Aborted path  %5d, step:%5d, t:%7.2e, h:%7.2e, |lambda(0)|: %7.2e\n',...
                                    ind, step, t, h, abs(y(n+1)))
                           end
                           abort_path = 1;
                        end
                    end
                end
                % check if path leads to an infinite eigenvalue
                if abort_inf>0 && abs(y(n+1))<abort_inf
                     % infinite eigenvalue
                     if disp_level==2
                        fprintf('Infinite path %5d, step:%5d, 1-t:%7.2e, h:%7.2e, |lambda(0)|: %7.2e\n',...
                            ind, step, 1-t, h, abs(y(n+1)))
                     end
                     converged(ind) = 2; 
                     abort_path = 1;
                end
            end
            tcell{ind} = tn;
            ycell{ind} = yn;
            X(:,ind) = yn(1:n,end);
            lambda(ind,:) = yn(n+1:end,end).';
            if t>=1-1e-15 && (converged(ind) == 0)
                converged(ind) = 1;
            end
            if disp_level==3
                fprintf('Fin path %5d, steps:%5d, minh: %7.2e,  1-t_end:%7.2e, |lam(0)|: %7.2e, status: %d \n',...
                    ind, step, minh(ind), 1-t, abs(y(n+1)), converged(ind))
            end
        end
    end  % parfoor loop

    lambdaN = lambda(:,2:end)./lambda(:,1);
    conv_ind = find(converged==1); % we do not repeat the computation for infinite eigenvalues
    inf_eigs = length(find(converged==2));

    % we repeat duplicate paths (details depend on the selected option)
    if strcmp(repeat_opt,'all')
        indices_dupl = conv_ind(find_duplicates(lambdaN(conv_ind,:),dist_tol)); 
    elseif strcmp(repeat_opt,'real')
        indices_dupl = conv_ind(rows_to_repeat(lambdaN(conv_ind,:),dist_tol,1e-8));
    else
        error('unknown option for repeat_eig')
    end
    
    % we also repeat paths that were aborted (and not detected as infinite eigenvalues)
    indices_fail = find(converged==0);
    fail = length(indices_fail);
    rep_indices = unique([indices_dupl indices_fail]);

    stat.follow_ind = follow_ind;
    stat.indices = rep_indices;
    stat.run = run;
    stat.converged = converged;

    if disp_level>0
        nconv = length(find(converged(follow_ind)==1)) + length(find(converged(follow_ind)==2));
        steps = cellfun(@numel, tcell);
        steps = steps(follow_ind);
        minstepused = min(minh(follow_ind));
        fprintf('Run %d, conv. %d/%d paths, fail %d, to repeat %d, inf %d; steps (avg %d, max %d, min %d); minh %6.1e, init %6.1e, max %6.1e; maxang %6.1e; innertol %6.1e, Newt %d\n',...
            run, nconv, length(follow_ind), fail, length(rep_indices), inf_eigs, round(mean(steps)),max(steps),min(steps),minstepused,stepsize,maxstepsize,maxangle,innertol,maxinnersteps)
    end
    if ~isempty(rep_indices)
        converged(rep_indices)=0;
        % we tighten criteria for the repeated curves
        stepsize = max(stepsize/100,minstepsize);
        maxstepsize = maxstepsize/3;
        innertol = max(innertol/sqrt(10),sqrt(10)*1e-15);
        maxinnersteps = max(maxinnersteps-1,4);
        maxangle = maxangle/3;
        maxsteps = 3*maxsteps;
    end
end

