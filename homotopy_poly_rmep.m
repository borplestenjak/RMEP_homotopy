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
%       - 'real' all real duplicates and complex duplicates without conjugate match
%   - use_pinv (0): solve Jacobian with pinv when the system is close to singular
%   - filter_res (1e-8): treshold for residual for finite eigenvalues
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

if isfield(opts,'stepsize'),      stepsize = opts.stepsize;           else,  stepsize = 1e-8;            end
if isfield(opts,'minstepsize'),   minstepsize = opts.minstepsize;     else,  minstepsize = 2*eps(class_t); end
if isfield(opts,'maxstepsize'),   maxstepsize = opts.maxstepsize;     else,  maxstepsize = 2.5e-2;       end
if isfield(opts,'maxsteps'),      maxsteps = opts.maxsteps;           else,  maxsteps = 2000;            end
if isfield(opts,'maxinnersteps'), maxinnersteps = opts.maxinnersteps; else,  maxinnersteps = 6;          end
if isfield(opts,'innertol'),      innertol = opts.innertol;           else,  innertol = 3e-13;           end
if isfield(opts,'finaltol'),      finaltol = opts.finaltol;           else,  finaltol = 1e-15;           end
if isfield(opts,'maxruns'),       maxruns = opts.maxruns;             else,  maxruns = 3;                end
if isfield(opts,'display'),       disp_level = opts.display;          else,  disp_level = 0;             end
if isfield(opts,'dist_tol'),      dist_tol = opts.dist_tol;           else,  dist_tol = 1e-6;            end
if isfield(opts,'fail_stepsize'), fail_stepsize = opts.fail_stepsize; else,  fail_stepsize = 1e-30;      end
if isfield(opts,'abort_inf'),     abort_inf = opts.abort_inf;         else,  abort_inf = 1e-15;          end
if isfield(opts,'maxangle'),      maxangle = opts.maxangle;           else,  maxangle = 1e-1;            end
if isfield(opts,'repeat_opt'),    repeat_opt = opts.repeat_opt;       else,  repeat_opt = 'all';         end
if isfield(opts,'use_pinv'),      use_pinv = opts.use_pinv;           else,  use_pinv = 0;               end
if isfield(opts,'filter_res'),    filter_res = opts.filter_res;       else,  filter_res = 1e2*sqrt(eps(class_t));  end

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
goal_eps = 10*eps(class_t); % how close to 1 must we get for convergence

deg = max(sum(suppA,2)); % degree of the RMEP
% support for homogeneous representations of RMEPS A(lambda) and B(lambda)
suppA_h = [deg-sum(suppA,2) suppA];
suppB_h = [deg-sum(suppB,2) suppB];

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
            Ma0 = zeros(m+2,n+kpar+1,class_t);
            b0 = zeros(m+2,1,class_t);
            Ca = zeros(m,n_cellA,class_t);   
            Cb = zeros(m,n_cellB,class_t);   
            C = cell(1,kpar+1);   
            % we follow homogeneous RMEP, we convert initial eigenvalue to homogeneous form
            lambda0 = [1; Lambda0(ind,:).'];
            lambda0 = lambda0/norm(lambda0);
            true_t = 0;
            t = 0;
            y = [X0(:,ind); lambda0]; % initial vector for the homotopy
            tn = t;
            yn = y;
            step = 0;
            h = stepsize;
            abort_path = 0;
            ort = 1; % orientation, we switch after true_t > 1/2
            ta = 0;
            tb = 1;
            facA = 1;
            facB = 0;
            facAp = 1;
            facBp = 1;
           
            % path following
            while (true_t<1-goal_eps) && (step<maxsteps) && (abort_path==0)

                if ort == 1 && true_t>0.5
                    t = 1-t;
                    ta = 1;
                    tb = 0;
                    ort = -1;
                end

                step = step + 1;
                tol = innertol;
                if ort==-1 && h >= t
                    h = t;
                    tol = finaltol; % we increase tolerance for the last step
                end

                facB = tb-ort*t;
                facA = ta+ort*t;
                if linearRMEP
                    W = 0;
                    Wt = 0;
                    for j = 1:kpar+1
                        C{j} = facB*B{j} + facA*A{j}; 
                        W = W + y(n+j)*C{j};
                        Wt = Wt + y(n+j)*(A{j} - B{j}); 
                    end
                    for j = 1:kpar+1
                        Ma0(1:m,n+j) = C{j}*y(1:n);
                    end
                else
                    Aval = eval_rmep(A,suppA_h,y(n+1:end),Z0);
                    Bval = eval_rmep(B,suppB_h,y(n+1:end),Z0);
                    W = facB*Bval + facA*Aval;  
                    Wt = Aval - Bval;
                
                    for j = 1:n_cellA
                        Ca(:,j) = A{j}*y(1:n);
                    end
                    for j = 1:n_cellB
                        Cb(:,j) = B{j}*y(1:n);
                    end
                    Ma0(1:m,n+1:end) = facB*eval_rmep_der(Cb,suppB_h,y(n+1:end),Z1b) + facA*eval_rmep_der(Ca,suppA_h,y(n+1:end),Z1a);
                end
                Ma0(1:m,1:n) = W; 
                Ma0(m+1,1:n) = 2*y(1:n)';
                Ma0(m+2,n+1:end) = 2*y(n+1:end)';
                b0(1:m) = -ort*Wt*y(1:n);
            
                % tangent vector (with t'=1) for the predictor
                [tang, rcnd] = linsolve(Ma0,b0); 
                if use_pinv && rcnd<1e-14
                    tang = pinv(Ma0)*b0;
                end
                noconv = 1;
                while noconv && (abort_path==0)
            
                    % compute predictor by Euler's method
                    yp = y + ort*h*tang;
                    tp = t + ort*h;
                    facBp = tb - ort*tp;  % 1 - true_t
                    facAp = ta + ort*tp;  % true_t

                    % initial residual
                    if linearRMEP
                        W = 0;
                        for j = 1:kpar+1
                            C{j} = facBp*B{j} + facAp*A{j}; 
                            W = W + yp(n+j)*C{j};
                        end
                    else
                        W = facBp*eval_rmep(B,suppB_h,yp(n+1:end),Z0) + facAp*eval_rmep(A,suppA_h,yp(n+1:end),Z0);
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
                            Ma(1:m,n+1:end) = facBp*eval_rmep_der(Cb,suppB_h,yp(n+1:end),Z1b) + facAp*eval_rmep_der(Ca,suppA_h,yp(n+1:end),Z1a);
                        end
                        Ma(m+1,1:n) = 2*yp(1:n)';
                        Ma(m+2,n+1:end) = 2*yp(n+1:end)';
                        % solve the system with the Jacobian matrix Ma
                        [dx, rcnd] = linsolve(Ma,-rhs);
                        if use_pinv && rcnd<1e-14
                            dx = -(pinv(Ma)*rhs);
                        end
                        yp = yp + dx; 
                        yp(1:n) = yp(1:n)/norm(yp(1:n));
                        yp(n+1:end) = yp(n+1:end)/norm(yp(n+1:end));
                        if linearRMEP
                            W = 0;
                            for j = 1:kpar+1
                                C{j} = facBp*B{j} + facAp*A{j}; 
                                W = W + yp(n+j)*C{j};
                            end
                        else    
                            W = facBp*eval_rmep(B,suppB_h,yp(n+1:end),Z0) + facAp*eval_rmep(A,suppA_h,yp(n+1:end),Z0);
                        end
                        rhs = [W*yp(1:n); yp(1:n)'*yp(1:n)-1; yp(n+1:end)'*yp(n+1:end)-1];
                    end

                    angle_vec = norm(y(1:n)-(y(1:n)'*yp(1:n))*yp(1:n));
                    angle_eig = norm(y(n+1:end)-(y(n+1:end)'*yp(n+1:end))*yp(n+1:end));
                    
                    if norm(rhs)<tol*normAB && angle_vec<maxangle && angle_eig<maxangle
                        % Newton method converged, angles are not too large, we take step forward
                        t = tp;
                        y = yp;
                        tn = [tn ta+ort*t];
                        yn = [yn y];
                        noconv = 0;
                        true_t = ta + ort*t;
                        % if convergence is very fast, we increase h 
                        if innerstep < maxinnersteps-1
                            h = min(2*h,maxstepsize);
                        end
                    else
                        % no convergence, we repeat predictor-corrector step with smaller h
                        h = h/2;
                        minh(ind) = min(minh(ind),h);
                        % if h is too small, we abort this path as we can not move forward in time anymore
                        if h<(t*minstepsize + fail_stepsize)
                           if disp_level>=2 
                                fprintf('Aborted path  %5d, step:%5d, ort: %3.1f, t:%7.2e, h:%7.2e, |lambda(0)|: %7.2e, rcnd: %7.2e\n',...
                                    ind, step, ort, t, h, abs(y(n+1)), rcnd)
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
                        fprintf('Infinite path %5d, step:%5d, ort: %3.1f, t:%7.2e, h:%7.2e, |lambda(0)|: %7.2e\n',...
                            ind, step, ort, t, h, abs(y(n+1)))
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
            if true_t>=1-goal_eps && (converged(ind) == 0) 
                converged(ind) = 1;
                % for all converged eigenvalues (homogeneous, could be finite or infinite) we compute 
                %  - the residual as an affine eigenvalue
                %  - the condition number as a homogeneous eigenvalue
                %  - the condition number as an affine eigenvalue
                lambdaA = yn(n+2:end,end).'/yn(n+1,end); % affine eigenvalue
                W = eval_rmep(A,suppA,lambdaA);
                restmp = norm(W*yn(1:n,end));
                [gm1,s1] = condeig_rmep(A,suppA,lambdaA);
                [gm2,s2] = condeig_rmep(A,suppA_h,yn(n+1:end,end).',false); % false because of homogeneous eigenvalues 
                rec_cond_all(ind,:) = [restmp gm1 s1 gm2 s2];
            end
            if disp_level>=4
                fprintf('Fin path %5d, steps:%5d, minh: %7.2e,  t_end:%7.2e, |lam(0)|: %7.2e, status: %d \n',...
                    ind, step, minh(ind), true_t, abs(y(n+1)), converged(ind))
            end
        end
    end  % parfoor loop

    lambdaA = lambda(:,2:end)./lambda(:,1); % affine eigenvalues
    conv_ind = find(converged==1); % indices of converged homogeneous eigenvalues

    rec_cond = rec_cond_all(conv_ind,:);

    % We divide converged paths into finite and infinite eigs
    % finite eigs are ones with small affine residual
    fin_ind = find(rec_cond(:,1)<filter_res);
    inf_ind = find(rec_cond(:,1)>=filter_res);
    fin_conv_ind = conv_ind(fin_ind); % indices of converged finite eigs
    inf_conv_ind = conv_ind(inf_ind); % indices of remaining converged eigs

    lambda_fin = lambdaA(fin_conv_ind,:); % converged finite affine eigs
    lambda_inf = lambda(inf_conv_ind,:);  % remaining converged eigs
    % we scale inf eigenvalues so that maximal element is 1 for easier search of duplicates
    for j = 1:length(inf_ind)
        [~,posmx] = max(abs(lambda_inf(j,2:end)));
        lambda_inf(j,:) = lambda_inf(j,:)/lambda_inf(j,posmx);
    end

    % we search for duplicates among finite eigenvalues
    fin_mult = find_duplicates(lambda_fin,dist_tol);       % indices of multiple finite eigenvalues
    fin_simple = setdiff(1:length(fin_conv_ind),fin_mult); % indices of simple finite eigenvalues
    % we repeat computation of all ill-conditioned simple eigenvalues
    fin_simple_ill_cond = fin_simple(find(rec_cond(fin_ind(fin_simple),3)>1/filter_res));  

    % in addition we repeat 
    %  - all multiple eigenvalues ('all')
    %  - just multiple eigenvalues that represent simple well-conditioned eigenvalues ('simple')
    %  - all multiple real eigenvalues and multiple complex without conjugate pairs ('real')
    if strcmp(repeat_opt,'all') 
        fin_repeat_dupl = fin_mult; 
    elseif strcmp(repeat_opt,'simple')
        fin_repeat_dupl = fin_mult(rec_cond(fin_ind(fin_mult),3)<1/filter_res);
    elseif strcmp(repeat_opt,'real')
        fin_repeat_dupl = rows_to_repeat(lambda_fin,dist_tol,filter_res);
    else
        error('unknown option for repeat_eig')
    end

    ind_repeat_fin = fin_conv_ind(union(fin_repeat_dupl,fin_simple_ill_cond));  % indices of finite eigenvalues to repeat

    % we search for duplicates among infinite eigenvalues
    inf_mult = find_duplicates(lambda_inf,dist_tol);
    inf_simple = setdiff(1:length(inf_conv_ind),inf_mult);
    % we repeat computation of all ill-conditioned infinite eigenvalues (in homogeneous setting)
    inf_simple_ill_cond = inf_simple(find(rec_cond(inf_ind(inf_simple),5)>1/filter_res));

    if strcmp(repeat_opt,'all') 
        inf_repeat_dupl = inf_mult; 
    elseif strcmp(repeat_opt,'simple')
        inf_repeat_dupl = inf_mult(rec_cond(inf_ind(inf_mult),3)<1/filter_res);
    elseif strcmp(repeat_opt,'real')
        inf_repeat_dupl = rows_to_repeat(lambda_inf,dist_tol,filter_res);
    else
        error('unknown option for repeat_eig')
    end
    
    ind_repeat_inf = inf_conv_ind(union(inf_repeat_dupl,inf_simple_ill_cond));  % indices of infinite eigenvalues to repeat
    
    indices_recomp = union(ind_repeat_fin,ind_repeat_inf); % indices of path to recompute

    % for some problems (arma, lti, ..) it is important to know the number of real finite solutions
    ind_real = find(vecnorm(double(imag(lambda_fin)),Inf,2)<filter_res);

    % we also repeat all aborted paths (and not detected as infinite eigenvalues)
    indices_fail = find(converged==0);
    fail = length(indices_fail);
    rep_indices = unique([indices_recomp indices_fail]);

    if nargout>4
        stat.converged = converged;     % status of convergence
        stat.rec_cond = rec_cond;       % residuals and condition number of all eigenvalues 
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
        minstepused = min(minh(follow_ind));
        inf_eigs = length(find(converged==2));
        if disp_level>1
           fprintf('Run %d',run)
        end
        fprintf(', conv %d, real %d, fin %d, to repeat %d, inf %d; steps (avg %d, max %d, min %d); minh %6.1e, h0 %6.1e, maxh %6.1e; maxfi %6.1e; intol %6.1e, Newt %d\n',...
            nconv, length(ind_real), length(fin_ind), length(rep_indices), inf_eigs, round(mean(steps)),max(steps),min(steps),minstepused, ...
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