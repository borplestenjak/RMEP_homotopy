function [DA,DB,indA,indB] = compare_sets(A,B,tol)

DA = [];
DB = [];
indA = [];
indB = [];

for j = 1:size(A,1)
    if ~is_in_set(B,A(j,:),tol)
        DA = [DA; A(j,:)];
        indA = [indA j];
    end
end

for j = 1:size(B,1)
    if ~is_in_set(A,B(j,:),tol)
        DB = [DB; B(j,:)];
        indB = [indB j];
    end
end

