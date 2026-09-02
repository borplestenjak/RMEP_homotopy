% Example four disc
% we compute H2 optimal model reducction 
% We get the same matrices as in the paper
% by Sibren Lagauw, Oscar Mauricio Agudelo, and Bart De Moor

b = [ 0.0448; 0.2368; 0.0013; 0.0211; 0.225; 0.0219 ];
a = [ -1.2024; 2.3675; -2.0039; 2.2337; -1.042; 0.8513];

addpath examples

red_order = 2;
[mat,supp] = h2sisored3_discrete(a,b,red_order);

opts = [];
opts.display = 1;
opts.maxruns = 3;
opts.abort_inf = 1e-4;  
opts.maxstepsize = 2e-1;
opts.maxangle = 1e-1;   
opts.stepsize = 1e-8;   
delta = 1e-5;

tstart = tic;
[pointsA,X] = poly_rect_multipareig_homotopy(mat,supp,opts);
tTrace = toc(tstart)
indh = find(vecnorm(imag(pointsA),Inf,2) < delta);
points = real(pointsA(indh, :))
X = X(:,indh);
for j=1:length(indh)
    X(:,j) = X(:,j)./X(1,j);
end
X
nsol2 = length(points)
disp(['Homotopy    solver required ', num2str(tTrace), 's'])

G_orig = tf(b',[1 a'],-1);
err = [];
ar = [];
br = [];
for j=1:length(indh)
    ar(j,:) = real([1 points(j,end:-1:1)]);
    br(j,:) = real(X(2:red_order+1,j).');
    G_red = tf(br(j,:),ar(j,:),-1);
    warning off
    err(j,1) = norm(G_orig-G_red,2);
    warning on
end

[tmperr,ord] = sort(err);
ind = find(tmperr<Inf);
sol_a = ar(ord(ind),:)
sol_b = br(ord(ind),:)
sol_err = err(ord(ind))
