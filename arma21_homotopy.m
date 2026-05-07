function [points,val,err,cand,simple] = arma21_homotopy(y,opts)

%ARMA21_HOMOTOPY returns stationary points for the ARMA(2,1) model
%
% [points, val, err, cand] = arma21(y,opts) returns real stationary points
% for the ARMA(2,1) model and objective function ||e||_2^2 
% 
% y(k) + alpha(1)*y(k-1) + alpha(2)*y(k-2) = e(k) + gamma(1)*e(k-1), k = 3,...,n
%
% for a given vector y of length n.
%
% Input:
%   - y : real vector of size n
%   - opts : options 
%
% Options in opts:
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - delta (sqrt(eps)): treshold for the extraction of real solutions
%   - all options for poly_rect_multipareig_homotopy
%
% Output:
%   - points : matrix with real critical points [alpha(1) gamma(1)]
%   - val: values of the objective function 
%   - err: minimal singular values used to verify the solution
%   - cand: matrix with all eigenvalues [alpha(1) gamma(1)]
%
% We find critical values as eigenvalues of a rectangular two-parameter
% eigenvalue problem (A + alpha(1)*B + alpha(2)*C + gamma(1)*D + gamma(1)^2*E)*z = 0,
% where A,B,C,D are matrices of size (4n-5)*(4n-7) and solve RMEP by the
% homotopy method

% Bor Plestenjak 2026

narginchk(1,2);

% Analyse user supplied options, if any
if nargin < 2, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(y);
end

if isfield(opts,'delta'),    delta   = opts.delta;    else,   delta = sqrt(eps(class_t));   end
if ~isfield(opts,'maxruns'),     opts.maxruns = 1;         end 
if ~isfield(opts,'maxstepsize'), opts.maxstepsize = 1e-1;  end 
if ~isfield(opts,'maxangle'),    opts.maxangle = 1e-1;     end 
if ~isfield(opts,'repeat_opt'),  opts.repeat_opt = 'real'; end 
if ~isfield(opts,'stepsize'),    opts.stepsize = 1e-12;    end 
if ~isfield(opts,'display'),     opts.display = 2;         end 

y = y(:);
if ~isa(y,class_t), y = numeric_t(y,class_t); end

% we build (4n-5)x(4n-7) matrices such that stationary points are eigenvalues 
% of the rectangular MEP (A + alfa1*B + alfa2*C + gamma1*D + gamma1^2*E) 
[M000,M100,M010,M001,M002] = ARMA21_mat(y);
A = {M000,M100,M010,M001,M002};
suppA = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 2];
[m,n] = size(M000);

[lambda,X,lambdaT,XT] = poly_rect_multipareig_homotopy(A,suppA,opts);
neig = size(lambda,1);

alpha1 = lambda(:,1);
alpha2 = lambda(:,2);
gamma = lambda(:,3);

msvd = zeros(neig,1); 
Y0 = zeros(m,n,class_t);
for k = 1:neig
    W = eval_rmep(A,suppA,lambda(k,:),Y0);
    msvd(k,:) = min(svd(W));
end

ind = find(vecnorm(imag(lambda),Inf,2)<delta);
alpha1RR = alpha1(ind);
alpha2RR = alpha2(ind);
gammaRR = gamma(ind);
sigmaRR = [];
errRR = msvd(ind);
for k = 1:length(gammaRR)
    sigmaRR(k,1) = arma21_err(y,real(alpha1RR(k)),real(alpha2RR(k)),real(gammaRR(k)),class_t);
end

points = real([alpha1RR alpha2RR gammaRR]);
val = sigmaRR;
cand = [alpha1 alpha2 gamma];
err = errRR;

% find simple finite solutions
mult = []; 
cn = []; 
for j=1:neig
    [gm,S] = condeig_rmep(A,suppA,lambda(j,:)); 
    mult(j,1) = gm; 
    cn(j,1) = S; 
end
ind_simple = find(mult==1);
simple = lambda(ind_simple,:);

end

function e = arma21_err(y,alfa1,alfa2,gamma,class_t)

    N = length(y);
    TC = diag(gamma*ones(N-2,1,class_t))+diag(ones(N-3,1,class_t),1); TC(N-2,N-1) = 1;
    TA = [diag(alfa2*ones(N-2,1,class_t)) zeros(N-2,2,class_t)]+[zeros(N-2,1,class_t) diag(alfa1*ones(N-2,1,class_t)) zeros(N-2,1,class_t)] + [zeros(N-2,2,class_t) diag(ones(N-2,1,class_t))];
    err = pinv(TC)*TA*y;
    e = norm(err)^2;

end

