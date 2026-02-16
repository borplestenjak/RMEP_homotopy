%% Define problem setup:
m = 4;
l = 5;
k = l + m - 1;

%% Build coefficient matrices:

M = cell(m + 1, 1);
A = [diag(randn(l,1)) + 1i*diag(randn(l,1)); zeros(m - 1,l)] + ...
    [zeros(m - 1, l); diag(randn(l,1)) + 1i*diag(randn(l, 1))];
M{1} = A; % make sure that A has one or two nonzero entries per row!
disp(A)
for i = 1:m
    M{i + 1} = -[zeros(i - 1,l); diag(ones(l, 1)); zeros(m - i,l)];
end

%% Test both approaches:

% Run default approach:
tic
disp("Default approach:")
W = cell(m,m + 1);
for i = 1:m
    W{i,1} = M{1};
    for j = 2:m+1
        W{i,j} = -M{j};
    end
end
TR = right_compression_matrix(l,m);
TL = left_compression_matrix(k,m);
DeltaTMP = multipar_delta(W);
Delta = cell(1, m + 1);
for i = 1:m + 1
    Delta{i} = TL*DeltaTMP{i}*TR;
end
t1 = toc;
disp(t1)

% Run new approach:
disp("New approach:")
tic
D = startdelta(M{1});
t2 = toc;
disp(t2)
