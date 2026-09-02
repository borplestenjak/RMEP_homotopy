b = [-2.9239 -39.5525 -97.5270 -147.1508];
a = [11.9584 43.9119 73.6759 44.3821];

a = a(:);
b = b(:);

addpath examples

red_order = 2;
[mat,supp] = h2sisored3(a,b,red_order);

scale_gamma = sqrt(norm(mat{1})/norm(mat{6}));
scale_delta = 2/(norm(mat{1})+norm(mat{2})*scale_gamma);
mats = mat;
mats{1} = scale_delta*mat{1};
for j = 2:3
    mats{j} = scale_gamma*scale_delta*mat{j};
end
for j = 4:6
    mats{j} = scale_gamma^2*scale_delta*mat{j};
end  

opts = [];
opts.display = 1;
opts.maxruns = 3;
opts.abort_inf = 1e-4;  
opts.maxstepsize = 2e-1;
opts.maxangle = 1e-1;   
opts.stepsize = 1e-8;   
delta = 1e-5;

tstart = tic;
[points,X] = poly_rect_multipareig_homotopy(mats,supp,opts);
points = points*scale_gamma;
tTrace = toc(tstart)
indh = find(vecnorm(imag(points),Inf,2) < delta);
points = real(points(indh, :))
X = X(:,indh);
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
for j=1:length(indh)
    ar(j,:) = real([1 points(j,end:-1:1)]);
    br(j,:) = real(X(2:red_order+1,j).');
    G_red = tf(br(j,:),ar(j,:));
    warning off
    err(j,1) = norm(G_orig-G_red,2)/normG;
    warning on
end

[tmperr,ord] = sort(err);
ind = find(tmperr<Inf);
sol_a = ar(ord(ind),:)
sol_b = br(ord(ind),:)
sol_err = err(ord(ind))




