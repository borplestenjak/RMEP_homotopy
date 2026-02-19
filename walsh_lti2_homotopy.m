function [points,val,err,cand,simple] = walsh_lti2_homotopy(y,opts)

%WALSH_LTI2_HOMOTOPY returns stationary points for the LTI(2) model
%
% [points, val, err, cand] = walsh_lti2_homotopy(y,opts) returns real stationary points
% for the LTI(2) model and objective function ||e||_2^2, where we are
% looking for parameters alpha(1) and alpha(2) for the best 2-norm
% approximation of y by z, such that its elements satisfy the difference
% equation 
%
% z(k+2) + alpha(1)*z(k+1) + alpha(2)*z(k) = 0, k = 1,...,n-2
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
% cubic eigenvalue problem and solve RMEP by the homotopy method

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
if ~isfield(opts,'repeat_opt'),  opts.repeat_opt = 'all'; end 
if ~isfield(opts,'stepsize'),    opts.stepsize = 1e-12;    end 
if ~isfield(opts,'display'),     opts.display = 1;         end 

y = y(:);
if ~isa(y,class_t), y = numeric_t(y,class_t); end

% we build matrices such that stationary points are eigenvalues 
% of the rectangular MEP 
[A,suppA] = lsrwalsh(y,2);
[m,n] = size(A{1});

[lambda,X,lambdaT,XT] = poly_rect_multipareig_homotopy(A,suppA,opts);
neig = size(lambda,1);

alpha1 = lambda(:,1);
alpha2 = lambda(:,2);

msvd = zeros(neig,1); 
Y0 = zeros(m,n,class_t);
for k = 1:neig
    W = eval_rmep(A,suppA,lambda(k,:),Y0);
    msvd(k,:) = min(svd(W));
end

ind = find(vecnorm(imag(lambda),Inf,2)<delta);
alpha1RR = alpha1(ind);
alpha2RR = alpha2(ind);
msvdRR = msvd(ind);
sigmaRR = zeros(length(alpha2RR),1);
for k = 1:length(alpha2RR)
    sigmaRR(k,1) = LTI2_err(y,real(alpha1RR(k)),real(alpha2RR(k)),class_t);
end

points = real([alpha1RR alpha2RR]);
val = sigmaRR;
cand = [alpha1 alpha2];
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

function e = LTI2_err(y,alpha1,alpha2,class_t)

    N = length(y);

    tmp1 = [alpha2*eye(N-2,class_t) zeros(N-2,2,class_t)];
    tmp2 = [zeros(N-2,1,class_t) alpha1*eye(N-2,class_t) zeros(N-2,1,class_t)];
    tmp3 = [zeros(N-2,2,class_t) eye(N-2,class_t)];
    TA = tmp1 + tmp2 + tmp3;
    DC = TA*TA';

    err = TA'*(DC\(TA*y));
    e = norm(err)^2;

end

