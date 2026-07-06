function [lambda,X,tcell,ycell,stat] = homotopy_poly_rmep(A,suppA,B,suppB,Lambda0,X0,opts)

% HOMOTOPY_POLY_RMEP   Solve a polynomial rectangular multiparameter 
% eigenvalue problem using homotopy
%
% [lambda,X,tn,yn,stat] = homotopy_poly_rmep(A,suppA,B,suppB,Lambda0,X0,opts) 
% returns homogeneous eigenvalues and eigenvectors of a polynomial rectangular 
% multiparameter eigenvalue problem (RMEP) 
%
% A(lambda) x = sum(j=1)^r1 (lambda.^suppA_r(j,:)*Aj) x = 0 
%
% using initial RMEP
%
% B(lambda) x = sum(j=1)^r2 (lambda.^suppB_r(j,:)*Bj) x = 0,
%
% where lambda=(1,lambda_1,...,lambda_k), with known initial eigenvalues 
% Lambda0 and eigenvectors X0
%
% Input:
%   - A,B: cell arrays containing coefficient matrices of a polynomial RMEP, 
%         each cell is an (n+k-1) x n matrix
%   - suppA,suppB: matrices with rows that contain degrees of monomials coefficient 
%         if suppA is empty, we assume that both RMEPs are linear
%   - Lambda0: rows are eigenvalues of the RMEP B (matrix of size m x k),
%         where m is the number of eigenvalues of B(lambda) x = 0
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
%   - maxruns (3): maximal number of attempts to follow curves
%   - display (0): display all iterations (>3),  infinite eigs (>2) aborted 
%        paths (>1), statistics of each run (1), no display (0) 
%   - dist_tol (1e-6): relative tolerance for checking if the solutions are distinct
%   - fail_stepsize (1e-30): abs. tolerance to abort path following if the stepsize becomes to small
%   - abort_inf (1e-15): tolerance on abs(lambda_0) to identify an infinite eigenvalue
%   - maxangle (1e-1): maximal angle (subspace) between eigenvectors and eigenvalues
%   - repeat_opt ('all'): how to select paths to repeat, 
%       - 'all': all duplicates, 
%       - 'simple': only duplicates that correspond to well-conditioned simple eigenvalues
%   - filter_res (1e-8): treshold for residual for finite eigenvalues
%   - gamma : random complex scalar in the homotopy (set to 1 to not use gamma)
%
% Output:
%   - lambda : matrix m x (k+1), each row is a homogeneous eigenvalue
%   - X : matrix n x m with right eigenvectors
%   - tcell: cell with all time steps 
%   - ycell: cell with all intermediate steps, first n elements represent eigenvector,
%         last k+1 elements represent eigenvalue 
%   - stat: several statistics

% Bor Plestenjak, 2026

% Validate number of input parameters
narginchk(6, 7);
if nargin < 7, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A{:});
end

if isfield(opts,'stepsize'),      stepsize = opts.stepsize;           else,  stepsize = numeric_t(1e-8,class_t); end
if isfield(opts,'minstepsize'),   minstepsize = opts.minstepsize;     else,  minstepsize = 2*eps(class_t); end
if isfield(opts,'maxstepsize'),   maxstepsize = opts.maxstepsize;     else,  maxstepsize = 2.5e-2;       end
if isfield(opts,'maxsteps'),      maxsteps = opts.maxsteps;           else,  maxsteps = 2000;            end
if isfield(opts,'maxinnersteps'), maxinnersteps = opts.maxinnersteps; else,  maxinnersteps = 6;          end
if isfield(opts,'innertol'),      innertol = opts.innertol;           else,  innertol = 1e3*eps(class_t);  end
if isfield(opts,'finaltol'),      finaltol = opts.finaltol;           else,  finaltol = 1e1*eps(class_t); end
if isfield(opts,'maxruns'),       maxruns = opts.maxruns;             else,  maxruns = 3;                end
if isfield(opts,'display'),       disp_level = opts.display;          else,  disp_level = 0;             end
if isfield(opts,'dist_tol'),      dist_tol = opts.dist_tol;           else,  dist_tol = 1e-6;            end
if isfield(opts,'fail_stepsize'), fail_stepsize = opts.fail_stepsize; else,  fail_stepsize = 1e2*eps(class_t)^2; end
if isfield(opts,'abort_inf'),     abort_inf = opts.abort_inf;         else,  abort_inf = 1e1*eps(class_t);  end
if isfield(opts,'maxangle'),      maxangle = opts.maxangle;           else,  maxangle = 1e-1;            end
if isfield(opts,'repeat_opt'),    repeat_opt = opts.repeat_opt;       else,  repeat_opt = 'all';         end
if isfield(opts,'filter_res'),    filter_res = opts.filter_res;       else,  filter_res = 1e2*sqrt(eps(class_t));  end
if isfield(opts,'goal_res'),      goal_eps = opts.goal_eps;           else,  goal_eps = 1e3*eps(class_t); end 
if isfield(opts,'gamma'),         gamma = opts.gamma;                 else,  gamma = exp(2*pi(class_t)*1i*rand(class_t)); end 

linearRMEP = isempty(suppA);
if linearRMEP 
    kpar = length(A) - 1;
else
    kpar = size(suppA,2);   % number of parameters
end
[m,n] = size(A{1});
neig = size(Lambda0,1); % number of paths to follow
tcell = cell(1,neig);
ycell = cell(1,neig);

jac_steps = zeros(1,neig); % number of systems that have to solved (predictor and correctors)
rep_indices = 1:neig; % indices of paths that we have to follow

n_cellA = numel(A);
n_cellB = numel(B);

% Make sure all inputs are of the same numeric type.
for k=1:numel(A)
    if ~isa(A{k}, class_t)
         A{k} = numeric_t(A{k},class_t);
    end
end
for k=1:numel(B)
    if ~isa(B{k}, class_t)
         B{k} = numeric_t(B{k},class_t);
    end
end
if ~isa(Lambda0,class_t), Lambda0 = numeric_t(Lambda0,class_t); end
if ~isa(X0,class_t), X0 = numeric_t(X0,class_t); end

stat = [];
stat.run = 0;

converged = zeros(1,neig);
minh = ones(1,neig,class_t);
run = 0;

deg = max(sum(suppA,2)); % degree of the RMEP
% support for homogeneous representations of RMEPS A(lambda) and B(lambda)
suppA_h = [deg-sum(suppA,2) suppA];
suppB_h = [deg-sum(suppB,2) suppB];

for k=1:numel(B) % introduce random scalar gamma
    B{k} = gamma*B{k};
end

tmpAB = [];
for k=1:numel(A)
    tmpAB = [tmpAB; norm(A{k})];
end
for k=1:numel(B)
    tmpAB = [tmpAB; norm(B{k})];
end
normAB = norm(tmpAB);
rec_cond_all = zeros(neig,5); % for residuals and condition numbers

while run<maxruns && length(find(converged==0))>0
    run = run + 1; 
    follow_ind = find(converged==0);
    minh(follow_ind) = 1;

    if disp_level>0
       fprintf('Run %d, %d paths', run, length(follow_ind))
       if disp_level>1
           fprintf('\n')
       end
    end

    stat.paths(run) = length(find(converged==0));
    
    % we follow paths in parallel, change parfor to for for debugging
    parfor ind = 1:neig
        if converged(ind)==0
            % initialization of data for path following
            Z0 = zeros(m,n,class_t);
            Z1a = zeros(n_cellA,kpar+1,class_t);
            Z1b = zeros(n_cellB,kpar+1,class_t);
            Ma = zeros(m+2,n+kpar+1,class_t);
            b0 = zeros(m+2,1,class_t);
            Ca = zeros(m,n_cellA,class_t);   
            Cb = zeros(m,n_cellB,class_t);   
            C = cell(1,kpar+1);   
            D = cell(1,kpar+1);   
            % we follow homogeneous RMEP, we convert initial eigenvalue to homogeneous form
            lambda0 = [1; Lambda0(ind,:).'];
            lambda0 = lambda0/norm(lambda0);
            t = numeric_t(0,class_t);
            y = [X0(:,ind); lambda0]; % initial vector for the homotopy
            tn = t;
            yn = y;
            step = 0;
            h = stepsize;
            abort_path = 0;
            jac_steps(ind) = 0;
            Aval = zeros(m,n,class_t);
            Bval = zeros(m,n,class_t);

            % initial Jacobian (t=0)   
            if linearRMEP
                W = numeric_t(0,class_t);
                Wt = numeric_t(0,class_t);
                for j = 1:kpar+1
                    C{j} = B{j}; 
                    D{j} = A{j}-B{j}; 
                    W = W + y(n+j)*B{j};
                end
            else
                Aval = eval_rmep(A,suppA_h,y(n+1:end),Z0,class_t);
                Bval = eval_rmep(B,suppB_h,y(n+1:end),Z0,class_t);
                W = Bval; 
            end

            % path following
            while (t<1-goal_eps) && (step<maxsteps) && (abort_path==0)

                step = step + 1;
                tol = innertol;
                if h > 1-t
                    h = 1-t;
                    tol = finaltol; % we increase tolerance for the last step
                end

                % Jacobian for the tangent vector
                Ma(1:m,1:n) = W;
                if linearRMEP
                    Wt = numeric_t(0,class_t);
                    for j = 1:kpar+1
                         Ma(1:m,n+j) = C{j}*y(1:n);
                        Wt = Wt + y(n+j)*D{j}; 
                    end
                else
                    for j = 1:n_cellA
                        Ca(:,j) = A{j}*y(1:n);
                    end
                    for j = 1:n_cellB
                        Cb(:,j) = B{j}*y(1:n);
                    end
                    Ma(1:m,n+1:end) = (1-t)*eval_rmep_der(Cb,suppB_h,y(n+1:end),Z1b,class_t) + t*eval_rmep_der(Ca,suppA_h,y(n+1:end),Z1a,class_t);
                    Wt = Aval - Bval;                    
                end
                Ma(m+1,1:n) = 2*y(1:n)';
                Ma(m+2,n+1:end) = 2*y(n+1:end)';
                b0(1:m) = -Wt*y(1:n);

                % tangent vector (with t'=1) for the predictor
                jac_steps(ind) = jac_steps(ind) + 1;
                [tang, rcnd] = linsolve(Ma,b0); 

                noconv = 1;
                while noconv && (abort_path==0)
            
                    % compute predictor by Euler's method
                    yp = y + h*tang;
                    tp = t + h;

                    % initial residual
                    if linearRMEP
                        W = numeric_t(0,class_t);
                        for j = 1:kpar+1
                            C{j} = (1-tp)*B{j} + tp*A{j}; 
                            W = W + yp(n+j)*C{j};
                        end
                    else
                        Aval = eval_rmep(A,suppA_h,yp(n+1:end),Z0,class_t);
                        Bval = eval_rmep(B,suppB_h,yp(n+1:end),Z0,class_t);
                        W = (1-tp)*Bval + tp*Aval; 
                    end
                    rhs = [W*yp(1:n); yp(1:n)'*yp(1:n)-1; yp(n+1:end)'*yp(n+1:end)-1];
                    innerstep = 0;

                    % compute corrector by Newton's method
                    while norm(rhs)>tol*normAB && innerstep<maxinnersteps
                        innerstep = innerstep + 1;
                        Ma(1:m,1:n) = W;
                        if linearRMEP
                            for j = 1:kpar+1
                                Ma(1:m,n+j) =  C{j}*yp(1:n);
                            end
                        else
                            for j = 1:n_cellA
                                Ca(:,j) = A{j}*yp(1:n);
                            end
                            for j = 1:n_cellB
                                Cb(:,j) = B{j}*yp(1:n);
                            end
                            Ma(1:m,n+1:end) = (1-tp)*eval_rmep_der(Cb,suppB_h,yp(n+1:end),Z1b,class_t) + tp*eval_rmep_der(Ca,suppA_h,yp(n+1:end),Z1a,class_t);
                        end
                        Ma(m+1,1:n) = 2*yp(1:n)';
                        Ma(m+2,n+1:end) = 2*yp(n+1:end)';
                        % solve the system with the Jacobian matrix Ma
                        jac_steps(ind) = jac_steps(ind) + 1;
                        [dx, rcnd] = linsolve(Ma,-rhs);
                        yp = yp + dx; 
                        yp(1:n) = yp(1:n)/norm(yp(1:n));
                        yp(n+1:end) = yp(n+1:end)/norm(yp(n+1:end));
                        % compute new residual H(t,yp)
                        if linearRMEP
                            W = numeric_t(0,class_t);
                            for j = 1:kpar+1
                                W = W + yp(n+j)*C{j};
                            end
                        else    
                            Aval = eval_rmep(A,suppA_h,yp(n+1:end),Z0,class_t);
                            Bval = eval_rmep(B,suppB_h,yp(n+1:end),Z0,class_t);
                            W = (1-tp)*Bval + tp*Aval; 
                        end
                        rhs = [W*yp(1:n); yp(1:n)'*yp(1:n)-1; yp(n+1:end)'*yp(n+1:end)-1];
                    end

                    angle_vec = norm(y(1:n)-(yp(1:n)'*y(1:n))*yp(1:n));
                    angle_eig = norm(y(n+1:end)-(yp(n+1:end)'*y(n+1:end))*yp(n+1:end));
                    
                    if norm(rhs)<tol*normAB && angle_vec<maxangle && angle_eig<maxangle
                        % Newton method converged, angles are not too large, we take step forward
                        t = tp;
                        y = yp;
                        tn = [tn t];
                        yn = [yn y];
                        noconv = 0;
                        % if convergence is very fast, we increase h 
                        if innerstep < maxinnersteps-1
                            h = min(sqrt(2)*h,maxstepsize);
                        end
                    else
                        % no convergence, we repeat predictor-corrector step with smaller h
                        h = h/2;
                        minh(ind) = min(minh(ind),h);
                        % if h is too small, we abort this path as we can not move forward in time anymore
                        if h<(t*minstepsize + fail_stepsize)
                           if disp_level>=2 
                                fprintf('Aborted path  %5d, step:%5d, t:%7.2e, h:%7.2e, |lambda(0)|: %7.2e, rcnd: %7.2e\n',...
                                    ind, step, t, h, abs(y(n+1)), rcnd)
                                fprintf('Parameters: %7.2e %7.2e %7.2e %7.2e %7.2e \n',norm(rhs),tol*normAB,angle_vec,angle_eig,maxangle)
                           end
                           abort_path = 1;
                        end
                    end
                end
                % check if path leads to an infinite eigenvalue
                if abort_inf>0 && abs(y(n+1))<abort_inf
                     % infinite eigenvalue
                     if disp_level>=3
                        fprintf('Infinite path %5d, step:%5d, t:%7.2e, h:%7.2e, |lambda(0)|: %7.2e\n',...
                            ind, step, t, h, abs(y(n+1)))
                     end
                     converged(ind) = 2;
                     abort_path = 1;
                end
            end
            % save results for the output
            tcell{ind} = tn;
            ycell{ind} = yn;
            X(:,ind) = yn(1:n,end);
            lambda(ind,:) = yn(n+1:end,end).';
            if t>=1-goal_eps && (converged(ind) == 0) 
                converged(ind) = 1;
            else
                if disp_level>=2 && step==maxsteps
                   fprintf('Maxsteps path  %5d, step:%5d, t:%7.2e, h:%7.2e, |lambda(0)|: %7.2e, 1-t:%7.2e\n',...
                          ind, step, t, h, abs(y(n+1)), 1-t)
                end
            end    
            if disp_level>=4
                fprintf('Fin path %5d, steps:%5d, minh: %7.2e,  t_end:%7.2e, |lam(0)|: %7.2e, status: %d \n',...
                    ind, step, minh(ind), t, abs(y(n+1)), converged(ind))
            end
        end
    end  % parfoor loop

    % we divide converged homogeneous eigenvalues into finite and infinite
    conv_ind = find(converged==1); % indices of converged homogeneous eigenvalues
    fin_ind = find(abs(lambda(conv_ind,1))>filter_res);
    inf_ind = find(abs(lambda(conv_ind,1))<filter_res);
    % we normalize finite eigenvalues so that lambda_0 = 1 and infinite so
    % that maximal element is 1
    lambdaN = lambda(conv_ind,:);
    lambdaN(fin_ind,:) = lambdaN(fin_ind,:)./lambdaN(fin_ind,1);
    for j = 1:length(inf_ind)
        [~,posmx] = max(abs(lambdaN(inf_ind(j),2:end)));
        lambdaN(inf_ind(j),:) = lambdaN(inf_ind(j),:)/lambdaN(inf_ind(j),posmx);
    end
    % now we can find duplicate indices for new trial
    ind_mult = find_duplicates(lambdaN,dist_tol);    % indices of multiple eigenvalues


    if strcmp(repeat_opt,'simple')
        % we repat only multiple eigenvalues that represent simple well-conditioned eigenvalues ('simple')
        cnd = [];
        for j = 1:length(ind_mult)
            [gm2,s2] = condeig_rmep(A,suppA_h,lambda(conv_ind(ind_mult(j)),:).',false); % false because of homogeneous eigenvalues 
            cnd(j,1) = s2;
        end
        indices_recomp = conv_ind(cnd<1/filter_res);
    else    
        % we repeat all multiple eigenvalues
        indices_recomp = conv_ind(ind_mult);  % indices of eigenvalues to repeat
    end
        
    % for some problems (arma, lti, ..) it is important to know the number of real finite solutions
    ind_real = find(vecnorm(double(imag(lambdaN(fin_ind,:))),Inf,2)<filter_res);

    % we also repeat all aborted paths (and not detected as infinite eigenvalues)
    indices_fail = find(converged==0);
    fail = length(indices_fail);
    rep_indices = unique([indices_recomp(:); indices_fail(:)]);

    if nargout>4
        stat.converged = converged;     % status of convergence
        stat.rep_indices{run} = rep_indices; % indices to recompute in the next runs
        stat.maxstepsize(run) = maxstepsize;   
        stat.stepsize(run) = stepsize; %  
        stat.stepsize(run) = stepsize; %  
        stat.innerrol(run) = innertol;
        stat.maxangle(run) = maxangle;
        stat.maxsteps(run) = maxsteps;
        stat.maxinnersteps(run) = maxinnersteps;
        stat.run = run;
    end
    
    if disp_level>0
        nconv = length(find(converged(follow_ind)==1)) + length(find(converged(follow_ind)==2));
        steps = cellfun(@numel, tcell);
        steps = steps(follow_ind);
        jacsteps = jac_steps(follow_ind);
        minstepused = min(minh(follow_ind));
        inf_eigs = length(find(converged==2));
        if disp_level>1
           fprintf('Run %d',run)
        end
        fprintf(', conv %d, real %d, fin %d, to repeat %d, inf %d; steps %d (avg %d, max %d, min %d); jacsteps %d (avg %d, max %d, min %d) minh %6.1e, h0 %6.1e, maxh %6.1e; maxfi %6.1e; intol %6.1e, Newt %d\n',...
            nconv, length(ind_real), length(fin_ind), length(rep_indices), inf_eigs, sum(steps), round(mean(steps)),max(steps),min(steps), ...
            sum(jacsteps), round(mean(jacsteps)),max(jacsteps),min(jacsteps), minstepused, ...
            stepsize, maxstepsize, maxangle, innertol, maxinnersteps)
    end
    if ~isempty(rep_indices)
        converged(rep_indices)=0;
        % we tighten criteria for the repeated curves
        stepsize = max(stepsize/100,minstepsize);
        maxstepsize = maxstepsize/3;
        innertol = max(innertol/sqrt(10),goal_eps);
        maxinnersteps = max(maxinnersteps-1,4);
        maxangle = maxangle/3;
        maxsteps = 3*maxsteps;
    end
end


