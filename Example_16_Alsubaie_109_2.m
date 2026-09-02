% Continuous-time H2 model reduction example
% from Alsubaie's Phd, page 109, example F1
% n = 5, r = 2 and r = 3

b = [-1.2805 -6.2266 -12.8095 -9.3373];
a = [3.1855 8.9263 12.2936 3.1987];

a = a(:);
b = b(:);

addpath examples

for red_order=2:3

    fprintf('Reduction to order %d\n-----------------------\n',red_order)
    [mat,supp] = h2sisored3(a,b,red_order);
    
    opts = [];
    opts.display = 1;
    opts.maxruns = 3;
    opts.abort_inf = 1e-4;  
    opts.maxstepsize = 2e-1;
    opts.maxangle = 1e-1;   
    opts.stepsize = 1e-8;   
    delta = 1e-5;
    
    [points2,X2] = poly_rect_multipareig_homotopy(mat,supp,opts);
    points2 = points2; %*scale_gamma;
    tTrace = toc(tstart)
    indh = find(vecnorm(imag(points2),Inf,2) < delta);
    points2 = real(points2(indh, :))
    X2 = X2(:,indh);
    for j=1:length(indh)
        X2(:,j) = real(X2(:,j)./X2(1,j));
    end
    nsol2 = length(points2)
    disp(['Homotopy    solver required ', num2str(tTrace), 's'])
    X2
    
    G_orig = tf(b',[1 a']);
    normG = norm(G_orig);
    err = [];
    ar = [];
    br = [];
    for j=1:length(indh)
        ar(j,:) = real([1 points2(j,end:-1:1)]);
        br(j,:) = real(X2(2:red_order+1,j).');
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
    poles = [];
    for j=1:length(ind)
        poles(j,:) = roots(sol_a(j,:));
    end
    poles
end

