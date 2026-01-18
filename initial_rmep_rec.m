function [B,lambda,X] = initial_rmep_rec(n,k)

%INITIAL_RMEP_REC   Initial linear rectangular MEP for the homotopy
% [B,lambda,X] = initial_rmep_rec(n,k) returns matrices, eigenvalues 
% and eigenvectors of a RMEP, using recursive construction for four or more
% parameters
%
% (B0 + lambda(1) B1 + ... + lambda(k) Bk) x = 0 
%
% Input:
%   - n: number of columns in matrices B
%   - k: number of parameters
%
% Output:
%   - B : cell {B0,B1,...,Bk}
%   - lambda : matrix m x k, each row is an eigenvalue
%   - X : matrix n x m with right eigenvectors

a = randn(n,1)+1i*randn(n,1);
b = randn(n,1)+1i*randn(n,1);

B = cell(1,k+1);
B{1} = [diag(a);zeros(k-1,n)] + [zeros(k-1,n); diag(b)];
for j = 2:k-1
    d = randn(n,1)+1i*randn(n,1);
    B{1} = B{1} + [zeros(j-1,n); diag(d); zeros(k-j,n)];
end
for j = 1:k
    B{j+1} = [zeros(j-1,n); eye(n); zeros(k-j,n)];
end

nsol = nchoosek(n+k-1,k);
lambda = zeros(nsol,k);
X = zeros(n,nsol);
ind_col = 0;
z0 = zeros(n,1);
A = cell(1,k-1);

for i = 1:n
    z = z0;
    z(i) = 1;
    for j = i:n
        lm = -a(i);
        um = -b(j);
        A{1} = B{1}(i+1:j+k-2,i:j) + lm*B{2}(i+1:j+k-2,i:j) + um*B{k+1}(i+1:j+k-2,i:j);
        for r = 2:k-1
            A{r} = B{r+1}(i+1:j+k-2,i:j);
        end
        m = j-i+1;
        kk = k-2;
        if i==j
            eta = -A{1}.';
            Z = 1;
        else
            choice = best_rec_solver(m,kk);
            if choice==1
                [eta,Z] = rect_multipareig(A);
            elseif choice==2
                [eta,Z] = rect_multipareig_macaulay(A);
            else
                [eta,Z] = rect_multipareig_homotopy(A);
            end
        end
        for r = 1:size(eta,1)
            ind_col = ind_col + 1;
            z = z0;
            z(i:j) = Z(:,r);
            lambda(ind_col, :)= [lm eta(r,:) um];
            X(:,ind_col) = z;
        end
    end
end

function choice = best_rec_solver(m,k)

% based on m and k we select the most efficient solver for the 
% subproblem RMEP 

% 1: rect_multipareig
% 2: rect_multipareig_Macaulay
% 3: rect_multipareig_homotopy

choice = 3; % default choice
if k==2
    choice = 1;
end
if k==3 && m<40
    choice = 1;
end
if k==4 && m<17
    choice = 1;
end
if (k==5 && m<6)
    choice = 2;
end
if (k==6 && m<=4) || (k==7 && m<=4)
    choice = 2;
end
if (k>=8 && m<=3)
    choice = 2;
end

