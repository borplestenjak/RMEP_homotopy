function [DA,DB] = compare_sets(A,B,tol)

DA = [];
DB = [];

for j = 1:size(A,1)
    if ~is_in_set(B,A(j,:),tol)
        DA = [DA; A(j,:)];
    end
end

for j = 1:size(B,1)
    if ~is_in_set(A,B(j,:),tol)
        DB = [DB; B(j,:)];
    end
end

