% Model reduction example FOM-2 from 
% Gugercin, Antoulas, Beaties
% n = 7, r = 3


b = [2 11.5 57.75 178.625 345.5 323.625 94.5];
a = [10 46 130 239 280 194 60];

a = a(:);
b = b(:);

addpath examples

red_order = 3;
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
[points2,X2] = poly_rect_multipareig_homotopy(mat,supp,opts);
tTrace = toc(tstart)
indh = find(vecnorm(imag(points2),Inf,2) < delta);
points = real(points2(indh, :))
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

