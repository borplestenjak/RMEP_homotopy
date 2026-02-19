% Example of a high multiplicity Van Vleck polynomial:
Q = cell(2, 1);
Q{1} = 0;
Q{2} = [0 0 0 0 1]; % Q2(z) = z^4;

[mat, supp] = heinestieltjes(Q, 3);
mep = mepstruct(mat, supp);

sol = macaulaylab(mep, verbose = true, clustering = false);

solMac = sol.num
solMPH = rect_multipareig_homotopy(mat)
solMP  = rect_multipareig(mat)