% This is model reduction Example 3 from Lagauw, Agudelo, De Moor
%
% n = 17 and m = 4
%
% S. Lagauw, O. M. Agudelo, and B. De Moor. Globally optimal SISO H2-norm
% model reduction using Walsh's theorem. IEEE Control Syst. Lett., 
% 7:1670–1675, 2023
    
load correctedh2example.mat
addpath examples

sys_ss = ss(A,B,C,D);
[sb,sa] = ss2tf(A,B,C,D); % transformation to transfer form

sb = sb(2:end);
n = length(sa)-1;
pden = n:-1:0;
pnum = length(sb)-1:-1:0;
gamma = 16; % scaling of coefficient matrices to reduce huge differences in size of the coefficients
sa2 = sa .* gamma.^pden;
sb2 = sb .* gamma.^pnum;

% Common normalization; does not change the transfer function
scale = max(abs(sa2(1)))
sa2 = sa2/scale;
sb2 = sb2/scale;

red_order = 2 % order of the reduced model (use 2,3, or 4)

sa2 = sa2(2:end);
sa2 = sa2(:);
sb2 = sb2(:);

G_orig = tf(sb2',[1 sa2']);

[mat,supp] = h2sisored3(sa2,sb2,red_order);

deltaopts = [];
opts.display = 1;
opts.maxruns = 2;
opts.maxstepsize = 2.5e-2;
opts.maxangle = 5e-2;   
opts.abort_inf = 1e-8; 
opts.stepsize = 1e-5;
opts.filter_res = 1e-5;
delta = 1e-5;
rng(1)

tstart = tic;
[points,X2,lambdaT,XT] = poly_rect_multipareig_homotopy(mat,supp,opts);
tTrace = toc(tstart);
indh = find(vecnorm(imag(points),Inf,2) < delta);
points1 = real(points(indh, :));
nsol1 = length(points1)
X = X2(:,indh);
for j=1:length(indh)
    X(:,j) = real(X(:,j)./X(1,j));
end
pointsR = points1;
disp(['Homotopy    solver required ', num2str(tTrace), 's'])

normG = norm(G_orig,2);
err = [];
ort = [];
ar = [];
br = [];
for j = 1:length(indh)
    ar(j,:) = real([1 pointsR(j,end:-1:1)]);
    br(j,:) = real(X(2:red_order+1,j).');
    G_red = tf(br(j,:),ar(j,:));
    warning off
    normRed = norm(G_red,2);
    normDif = norm(G_orig-G_red,2);
    err(j,1) = normDif/normG;
    ort(j,1) = normDif^2-(normG^2-normRed^2); % Should be zero for a stationary point
    warning on
end

[tmperr,ord] = sort(err);
ind = find(tmperr<Inf);
sol_a = ar(ord(ind),:)
sol_b = br(ord(ind),:)
sol_err = err(ord(ind));
sol_ort = ort(ord(ind));
sol_err_ort = [sol_err sol_err*normG*sqrt(gamma) sol_ort]
poles = [];
for j=1:length(ind)
    poles(j,:) = roots(sol_a(j,:));
end
poles = poles*gamma
