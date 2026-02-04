function [points,val,err,cand,simple] = arma11_homotopy(y,opts)

%ARMA11_HOMOTOPY returns stationary points for the ARMA(1,1) model
%
% [points, val, err, cand, simple] = arma11_homotopy(y,opts) returns real 
% stationary points for the ARMA(1,1) model and objective function ||e||_2^2 
% 
% y(k) + alpha(1)*y(k-1) = e(k) + gamma(1)*e(k-1), k = 2,...,n
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
%   - simple: matrix with simple eigenvalues 
%
% We find critical values as eigenvalues of a rectangular two-parameter
% eigenvalue problem (A + alpha(1)*B + gamma(1)*C + gamma(1)^2*D)*z = 0,
% where A,B,C,D are matrices of size (3n-1)*(3n-2) and solve RMEP by the
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
if ~isfield(opts,'display'),     opts.display = 1;         end 

y = y(:);
if ~isa(y,class_t), y = numeric_t(y,class_t); end

% we build (3n-1)x(3n-2) matrices such that stationary points are eigenvalues 
% of the rectangular MEP (M00 + lambda1*M10 + lambda2*M01 + lambda2^2*M02) 
[M00,M10,M01,M02] = ARMA11_mat(y);
A = {M00,M10,M01,M02};
suppA = [0 0; 1 0; 0 1; 0 2];
[m,n] = size(M00);

[lambda,X,lambdaT,XT] = poly_rect_multipareig_homotopy(A,suppA,opts);
neig = size(lambda,1);

alpha = lambda(:,1);
gamma = lambda(:,2);

msvd = zeros(neig,1); 
Y0 = zeros(m,n,class_t);
for k = 1:neig
    W = eval_rmep(A,suppA,lambda(k,:),Y0);
    msvd(k,:) = min(svd(W));
end

ind = find(vecnorm(imag(lambda),Inf,2)<delta);
alphaRR = alpha(ind);
gammaRR = gamma(ind);
msvdRR = msvd(ind);
sigmaRR = zeros(length(gammaRR),1);
for k = 1:length(gammaRR)
    sigmaRR(k,1) = arma11_err(y,real(alphaRR(k)),real(gammaRR(k)),class_t);
end

points = real([alphaRR gammaRR]);
val = sigmaRR;
cand = [alpha gamma];
err = msvdRR;

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

function e = arma11_err(y,alpha,gamma,class_t)
% returns value of the objective function for the ARMA(1,1) model

    N = length(y);
    TC = diag(gamma*ones(N-1,1,class_t))+diag(ones(N-2,1,class_t),1); TC(N-1,N) = 1;
    TA = diag(alpha*ones(N-1,1,class_t))+diag(ones(N-2,1,class_t),1); TA(N-1,N) = 1;
    err = pinv(TC)*TA*y;
    e = norm(err)^2;

end


