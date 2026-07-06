function [lambda,X,stat] = homotopy_poly_rmep_gpu(A,suppA,B,suppB,Lambda0,X0,opts)

% HOMOTOPY_RMEP_GPU  Solve a rectangular multiparameter 
% eigenvalue problem using homotopy
%
% [lambda,X,stat] = homotopy_rmep_gpu(A,suppA,B,suppB,Lambda0,X0,opts) 
% returns homogeneous eigenvalues and eigenvectors of a rectangular 
% multiparameter eigenvalue problem (RMEP) 
%
% (A{1} + lambda(1) A{2} x + ... + lambda(k) A{k+1}) x = 0 
%
% using initial RMEP
%
% (B{1} + lambda(1) B{2} x + ... + lambda(k) B{k+1}) x = 0 
%
% where lambda=(1,lambda_1,...,lambda_k), with known initial eigenvalues 
% Lambda0 and eigenvectors X0
%
% Input:
%   - A,B: cell arrays containing coefficient matrices of a RMEP, 
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

if isfield(opts,'stepsize'),      stepsize = opts.stepsize;           else,  stepsize = single(1e-5); end
if isfield(opts,'minstepsize'),   minstepsize = opts.minstepsize;     else,  minstepsize = 2*eps('single'); end
if isfield(opts,'maxstepsize'),   maxstepsize = opts.maxstepsize;     else,  maxstepsize = 1e-1;       end
if isfield(opts,'maxsteps'),      maxsteps = opts.maxsteps;           else,  maxsteps = 400;            end
if isfield(opts,'maxinnersteps'), maxinnersteps = opts.maxinnersteps; else,  maxinnersteps = 5;          end
if isfield(opts,'innertol'),      innertol = opts.innertol;           else,  innertol = 1e2*eps('single');  end
if isfield(opts,'finaltol'),      finaltol = opts.finaltol;           else,  finaltol = 1e1*eps('single'); end
if isfield(opts,'display'),       disp_level = opts.display;          else,  disp_level = 1;             end
if isfield(opts,'dist_tol'),      dist_tol = opts.dist_tol;           else,  dist_tol = 1e-6;            end
if isfield(opts,'fail_stepsize'), fail_stepsize = opts.fail_stepsize; else,  fail_stepsize = 1e2*eps('single')^2; end
if isfield(opts,'abort_inf'),     abort_inf = opts.abort_inf;         else,  abort_inf = 1e1*eps('single');  end
if isfield(opts,'maxangle'),      maxangle = opts.maxangle;           else,  maxangle = 2e-1;            end
if isfield(opts,'filter_res'),    filter_res = opts.filter_res;       else,  filter_res = 1e2*sqrt(eps('single'));  end
if isfield(opts,'goal_res'),      goal_eps = opts.goal_eps;           else,  goal_eps = 1e2*eps('single'); end % how close to 1 must we get for convergence
if isfield(opts,'minsol'),        minsol = opts.minsol;               else,  minsol = 20; end % how many eigenvalues can remain unconverged

gdev = gpuDevice;
if disp_level>0
    fprintf('Initial available memory: %.2f GB\n', gdev.AvailableMemory/1e9)
end

linearRMEP = isempty(suppA);
if linearRMEP 
    kpar = length(A) - 1;
else
    kpar = size(suppA,2);   % number of parameters
end
[m,n] = size(A{1});

neig = size(Lambda0,1); % number of paths to follow
rep_indices = 1:neig; % indices of paths that we have to follow

n_cellA = numel(A);
n_cellB = numel(B);

stat = [];
stat.run = 0;

converged = zeros(1,neig);
run = 0;

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

Ag = cellfun(@(M) gpuArray(single(M)), A, 'UniformOutput', false);
Bg = cellfun(@(M) gpuArray(single(M)), B, 'UniformOutput', false);
normABg = gpuArray(single(normAB));

lambda0_h = [ones(neig,1) single(Lambda0)];
lambda0_h = lambda0_h ./ vecnorm(lambda0_h,2,2); % homogenized initial eigenvalues

pgn = min(neig,2048); % nummber of pages in GPU for one batch

finished = false(1,neig,'gpuArray');          % path status flag
failed   = false(1,neig,'gpuArray');          % path status flag
h = gpuArray(stepsize*ones(1,neig,'single')); % stepsize for each path
tol = gpuArray(zeros(1,neig,'single'));       % tolerance for each path
y = gpuArray(zeros(n+kpar+1,neig,'single'));   % vector [x; lambda] for each path, size (n+k+1) x pgn
t = gpuArray(zeros(1,neig,'single'));          % position t for each path
Wt = zeros(m,neig,'single','gpuArray');

newton_steps = zeros(1,pgn,'uint8','gpuArray');
rhs  = zeros(m+2,1,pgn,'single','gpuArray');

y = single([X0; lambda0_h.']);              

active = ~finished & ~failed;
step = 0;

while step < maxsteps && gather(any(active)) 
    % in each step we move all nonconverged paths simultaneously on GPU    

    step = step + 1;
    active = ~finished & ~failed;
    tol(t+h>1) = finaltol;
    h_eff = min(h,1-t) .* single(active);

    tol(:) = innertol;

    for batch = 1:pgn:neig
        pA = batch;
        pB = min(batch+pgn-1,neig);
        indAB = pA:pB;
        npages = pB-pA+1;
        %if disp_level>0
        %   fprintf('Step %d, batch %d to %d / %d \n', step, pA,pB,neig)
        %end
        

        % 1) compute tangent vectors for all paths
        % 1a) compute matrices for linear systems
        % -------------------------------------------------------------------
        % Ma  = zeros(m+2,n+kpar+1,neig,'single','gpuArray');
        if linearRMEP 
            Ma0 = page_Jacobian(Ag,Bg,t(pA:pB),y(:,pA:pB),m,n,kpar,npages);
            Wt(:,1:npages) = 0;
            for j = 1:kpar+1
                Wt(:,1:npages) = Wt(:,1:npages) + y(n+j,pA:pB) .* ((Ag{j}-Bg{j})*y(1:n,pA:pB));
            end
            rhs  = zeros(m+2,1,npages,'single','gpuArray');
            rhs(1:m,1,1:npages) = reshape(-Wt(:,1:npages),m,1,npages);
        else
            [Ma0,~,rhs] = page_Jacobian_poly(Ag,suppA_h,Bg,suppB_h,t(pA:pB),y(:,pA:pB),m,n,kpar,npages);
        end
    
        % 1b) solve linear systems
        % -------------------------------------------------------------------
        tang = pagefun(@mldivide,Ma0,rhs);
    
        % 1c) predictor values
        % -------------------------------------------------------------------
        yp = y(:,pA:pB) + reshape(tang,n+kpar+1,npages).*h_eff(pA:pB);
        tp = t(pA:pB) + h_eff(pA:pB);
    
        % 2) Newton corrector
        % initial residual
        if linearRMEP   
            [Ma,res] = page_Jacobian(Ag,Bg,tp(1:npages),yp(:,1:npages),m,n,kpar,npages);
        else
            [Ma,res] = page_Jacobian_poly(Ag,suppA_h,Bg,suppB_h,tp(1:npages),yp(:,1:npages),m,n,kpar,npages);
        end
        newton_steps(1:npages) = uint8(0);

        normres = reshape(pagefun(@norm,res),1,npages); 
        active_Newton = normres > tol(indAB)*normABg;
        inner = 1;
        while inner < maxinnersteps && gather(any(active_Newton))   
            newton_steps(active_Newton) = newton_steps(active_Newton) + uint8(1);
            dx = pagefun(@mldivide,Ma,-res);
            mask = single(active_Newton);
            yp = yp + reshape(dx,n+kpar+1,npages) .* mask; 
            yp(1:n,active_Newton) = yp(1:n,active_Newton) ./ sqrt(sum(abs(yp(1:n,active_Newton)).^2,1));
            yp(n+1:end,active_Newton) = yp(n+1:end,active_Newton) ./ sqrt(sum(abs(yp(n+1:end,active_Newton)).^2,1));
            if linearRMEP   
                [Ma,res] = page_Jacobian(Ag,Bg,tp,yp,m,n,kpar,npages);
            else
                [Ma,res] = page_Jacobian_poly(Ag,suppA_h,Bg,suppB_h,tp,yp,m,n,kpar,npages);
            end
            inner = inner + 1;
            normres = reshape(sqrt(sum(abs(res).^2,1)),1,npages); 
            active_Newton = normres > tol(indAB)*normABg;
        end
        % check angle
        dotprod_vec = sum(conj(yp(1:n,:)).*y(1:n,pA:pB),1);      % 1 x p
        angle_vec = sqrt(sum(abs(y(1:n,pA:pB) - yp(1:n,:).*dotprod_vec).^2,1));
        dotprod_eig = sum(conj(yp(n+1:end,:)).*y(n+1:end,pA:pB),1);      % 1 x p
        angle_eig = sqrt(sum(abs(y(n+1:end,pA:pB) - yp(n+1:end,:).*dotprod_vec).^2,1));
        active_Newton = active_Newton | (angle_vec > maxangle) | (angle_eig > maxangle);

        % update stepsizes for next step
        accept_ind = active(indAB) & ~active_Newton;
        y(:,indAB(accept_ind)) = yp(:,accept_ind);
        t(indAB(accept_ind)) = tp(accept_ind);
        enlarge_step = accept_ind & (newton_steps(1:npages)<2);
        h(indAB(enlarge_step)) = min(1.25*h(indAB(enlarge_step)),maxstepsize);
        finished = finished | (t >= 1 - goal_eps);
        
        % reject = active & ~accept;
        h(indAB(active_Newton)) = h(indAB(active_Newton))/2;
        failed = failed | (h < fail_stepsize);
    end
        
        active = ~finished & ~failed;

        if disp_level>0
            fprintf('Newton finished step %d, finished: %d failed: %d  active %d, max_t %7.3e, min_t %7.3e Mem: %.2f GB, \n',step,gather(sum(finished)),gather(sum(failed)),gather(sum(active)),gather(max(t)),gather(min(t)), gdev.AvailableMemory/1e9)
        end
        if gather(sum(active))<minsol
            maxsteps = min(maxsteps,step+5);
            minsol = 0;
        end
    
end
    
lambda = gather(y(n+1:end,:)).';
X = gather(y(1:n,:));

end

% Helper functions

function Y = page_times_vec(W,X)
 % multiplication of page_matrix with page_vector
 %      W of size m x n x p
 %      X of size n x p
 % returns an m x p matrix with columns W(:,:,i)*X(,:) for i=1,..,p
    [m,n,p] = size(W);
    X3 = reshape(X,n,1,p);
    Y3 = pagefun(@mtimes,W,X3);
    Y = reshape(Y3,m,p);
end

function [Jac, res, rhs] = page_Jacobian(Ag,Bg,t,y,m,n,kpar,pgn)
% page Jacobian matrices and residuals for Newton's method 
    Jac = zeros(m+2,n+kpar+1,pgn,'single','gpuArray');
    % res  = zeros(m+2,1,pgn,'single','gpuArray');
    W = zeros(m,n,pgn,'single','gpuArray');
    D = zeros(m,kpar+1,pgn,'single','gpuArray');
    for j = 1:kpar+1
        Cj = reshape((1-t),1,1,pgn).*Bg{j} + reshape(t,1,1,pgn).*Ag{j};
        W = W + reshape(y(n+j,:),1,1,pgn).*Cj;
        CjX = ((1-t).*(Bg{j}*y(1:n,:))) +(t.*(Ag{j}*y(1:n,:)));
        D(:,j,:) = reshape(CjX,m,1,pgn);
    end
    Jac(1:m,1:n,:) = W;
    Jac(1:m,n+1:end,:) = D;
    Jac(m+1,1:n,:) = conj(reshape(2*y(1:n,:),1,n,pgn));
    Jac(m+2,n+1:end,:) = conj(reshape(2*y(n+1:end,:),1,kpar+1,pgn)); 
    res(1:m,1,:) = page_times_vec(W,y(1:n,:));
    res(m+1,1,:) = reshape(sum(abs(y(1:n,:)).^2,1)-1,1,1,pgn);
    res(m+2,1,:) = reshape(sum(abs(y(n+1:end,:)).^2,1)-1,1,1,pgn);   
end

function [Jac, res, rhs] = page_Jacobian_poly(Ag,suppA_h,Bg,suppB_h,t,y,m,n,kpar,pgn)
% page Jacobian matrices and residuals for Newton's method 
    Jac = zeros(m+2,n+kpar+1,pgn,'single','gpuArray');
    rhs  = zeros(m+2,1,pgn,'single','gpuArray');
    [WA,DA] = eval_poly_WD_batch(Ag,suppA_h,y,m,n,kpar,pgn);
    [WB,DB] = eval_poly_WD_batch(Bg,suppB_h,y,m,n,kpar,pgn);
    tt = reshape(t,1,1,pgn);
    W = tt.*WA + (1-tt).*WB;
    D = tt.*DA + (1-tt).*DB;
    Jac(1:m,1:n,:) = W;
    Jac(1:m,n+1:end,:) = D;
    Jac(m+1,1:n,:) = conj(reshape(2*y(1:n,:),1,n,pgn));
    Jac(m+2,n+1:end,:) = conj(reshape(2*y(n+1:end,:),1,kpar+1,pgn)); 
    res(1:m,1,:) = page_times_vec(W,y(1:n,:));
    res(m+1,1,:) = reshape(sum(abs(y(1:n,:)).^2,1)-1,1,1,pgn);
    res(m+2,1,:) = reshape(sum(abs(y(n+1:end,:)).^2,1)-1,1,1,pgn);   
    rhs(1:m,1,:) = page_times_vec(WB,y(1:n,:)) - page_times_vec(WA,y(1:n,:));
end 

function [W,D] = eval_poly_WD_batch(Cg,supp,y,m,n,kpar,pgn)

    r = numel(Cg);
    W = zeros(m,n,pgn,'single','gpuArray');
    D = zeros(m,kpar+1,pgn,'single','gpuArray');
    dmon = ones(1,pgn,'single','gpuArray');
    deg = max(max(supp));
    Pow = cell(kpar+1,deg+1);
    for j = 1:kpar+1
        Pow{j,1} = ones(1,pgn,'single','gpuArray');
        Pow{j,2} = y(n+j,:);
        for d = 2:deg 
            Pow{j,d+1} = Pow{j,d} .* y(n+j,:);
        end
    end

    for q = 1:r
        mon = Pow{1,supp(q,1)+1}; 
        for j = 2:kpar+1
            mon = mon .* Pow{j,supp(q,j)+1}; 
        end
        W = W + reshape(mon,1,1,pgn).*Cg{q};
        CX = Cg{q} * y(1:n,:);
        for j = 1:kpar+1
            a = supp(q,j);
            if a ~= 0
                dmon(:) = a; 
                for rr = 1:kpar+1
                    pow = supp(q,rr);
                    if rr == j
                        pow = pow - 1;
                    end
                    if pow ~= 0
                        dmon = dmon .* Pow{rr, pow+1}; 
                    end
                end
                D(:,j,:) = D(:,j,:) + reshape(dmon.*CX,m,1,pgn);
            end
        end
    end
end

function z = gpu_int_power(x,pow)

    if pow == 0
        z = ones(size(x),'single','gpuArray');
    elseif pow == 1
        z = x;
    elseif pow == 2
        z = x .* x;
    else
        z = ones(size(x),'single','gpuArray');
        for kk = 1:pow
            z = z .* x;
        end
    end

end