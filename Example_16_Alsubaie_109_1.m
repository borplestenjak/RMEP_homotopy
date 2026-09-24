% Continuous-time H2 model reduction example
% from Alsubaie's Phd, page 109, example F1
% n = 5, r = 2


b = [-2.9239 -39.5525 -97.5270 -147.1508];
a = [11.9584 43.9119 73.6759 44.3821];

a = a(:);
b = b(:);

addpath examples

red_order = 2;
[mat,supp] = h2sisored3(a,b,red_order);

opts = [];
opts.display = 1;
opts.maxruns = 3;
opts.abort_inf = 1e-4;  
opts.maxstepsize = 2e-1;
opts.maxangle = 1e-1;   
opts.stepsize = 1e-8;   
delta = 1e-5;

tstart = tic;
[points,X2] = poly_rect_multipareig_homotopy(mat,supp,opts);
points = points; 
tTrace = toc(tstart)
indh = find(vecnorm(imag(points),Inf,2) < delta);
points = real(points(indh, :))
X = X2(:,indh);
for j=1:length(indh)
    X(:,j) = real(X(:,j)./X(1,j));
end
nsol2 = length(points)
disp(['Homotopy    solver required ', num2str(tTrace), 's'])
X

G_orig = tf(b',[1 a']);
normG = norm(G_orig);
err = [];
ar = [];
br = [];
ort = [];
for j=1:length(indh)
    ar(j,:) = real([1 points(j,end:-1:1)]);
    br(j,:) = real(X(2:red_order+1,j).');
    G_red = tf(br(j,:),ar(j,:));
    warning off
    normRed = norm(G_red,2);
    normDif = norm(G_orig-G_red,2);
    err(j,1) = normDif/normG;
    ort(j,1) = normDif^2-(normG^2-normRed^2);
    warning on
end

[tmperr,ord] = sort(err);
ind = find(tmperr<Inf);
sol_a = ar(ord(ind),:)
sol_b = br(ord(ind),:)
sol_err = err(ord(ind));
sol_ort = ort(ord(ind));
sol_err_ort = [sol_err sol_ort]
poles = [];
for j=1:length(ind)
    poles(j,:) = roots(sol_a(j,:));
end
poles

