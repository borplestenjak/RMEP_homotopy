% Continuous-time H2 model reduction example
% from Alsubaie's Phd, page 110
% n = 5, r = 3

SA = [-12.9672  0.2462  7.6147 -18.1222  -8.1217;
       -0.6374 -2.8920  2.3592 -13.9644  -9.8803;
       -2.7476  0.3302  5.7949 -17.9580 -11.8268;
       -3.9880 -0.8316  8.0971 -23.7215 -14.8762;
        5.6756  0.8380 -7.3135  21.7049  12.1270];
Sb = [-0.3045; -0.3632;-0.0579; 1.0672; -0.6764];
Sc = [ 0.2525  -1.9107 -0.9154  0.9401  -1.7435];

sys_ss = ss(SA,Sb,Sc,0);
sys_tf = tf(sys_ss);
sb = sys_tf.Numerator{1};
sa = sys_tf.Denominator{1};
sb = sb(2:end);
sa = sa(2:end);

sa = sa(:);
sb = sb(:);

addpath examples

red_order = 3;
[mat,supp] = h2sisored3(sa,sb,red_order);

opts = [];
opts.display = 1;
opts.maxruns = 3;
opts.abort_inf = 1e-4;  
opts.maxstepsize = 2e-1;
opts.maxangle = 1e-1;   
opts.stepsize = 1e-8;   
opts.filter_res = 1e-6;
delta = 1e-5;

rng(1)
tstart = tic;
[points2,X2,lambdaT,XT] = poly_rect_multipareig_homotopy(mat,supp,opts);
points2 = points2; %*scale_gamma;
tTrace = toc(tstart)
indh = find(vecnorm(imag(points2),Inf,2) < delta);
points = real(points2(indh, :))
X2 = X2(:,indh);
for j=1:length(indh)
    X2(:,j) = real(X2(:,j)./X2(1,j));
end
nsol2 = length(points)
disp(['Homotopy    solver required ', num2str(tTrace), 's'])
X2

G_orig = tf(sb',[1 sa']);
normG = norm(G_orig);
err = [];
ar = [];
br = [];
for j=1:length(indh)
    ar(j,:) = real([1 points(j,end:-1:1)]);
    br(j,:) = real(X2(2:red_order+1,j).');
    G_red = tf(br(j,:),ar(j,:));
    warning off
    err(j,1) = norm(G_orig-G_red,2)/normG;
    warning on
end

[tmperr,ord] = sort(err);
ind = find(tmperr<Inf);ar
sol_a = ar(ord(ind),:)
sol_b = br(ord(ind),:)
sol_err = err(ord(ind))
poles = [];
for j=1:length(ind)
    poles(j,:) = roots(sol_a(j,:));
end
poles

