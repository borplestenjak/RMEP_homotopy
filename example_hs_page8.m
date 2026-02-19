% Example on page 8:
a = -2 - 3*1i;
b = -3 + 2*1i;
Q = cell(3, 1);
Q{1} = 0;
Q{2} = 0;
Q{3} = [a*b (a + b) (1 + a*b) (a + b) 1]; % Q3(z) = (z^2 + 1)(z + a)(z + b)

[mat, supp] = heinestieltjes(Q, 2);
mep = mepstruct(mat, supp);

sol = macaulaylab(mep, verbose = true, clustering = false);

solMac = sol.num
solMPH = rect_multipareig_homotopy(mat)
solMP  = rect_multipareig(mat)