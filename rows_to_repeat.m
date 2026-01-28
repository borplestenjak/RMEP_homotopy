function [ind,IA] = rows_to_repeat(A,tol,imagtol)

% ROWS_TO_REPEAT indices of rows to repeat in homotopy method

if nargin < 2 || isempty(tol),     tol = 1e-6; end
if nargin < 3 || isempty(imagtol), imagtol = 1e-10; end

% Group rows using uniquetol, 'ByRows' treats each row as a single point 
% in N-dimensional space, matrix must be real

M = [real(A) imag(A)];
rowNorms = sqrt(sum(abs(M).^2, 2));
rowNorms(rowNorms == 0) = 1;
M_norm = M ./ rowNorms; % normalized rows

[C1, ~, IC1] = uniquetol(M_norm, tol, 'ByRows', true, 'DataScale', 1);
[C2, ~, IC2] = uniquetol(rowNorms, tol, 'DataScale', 1);

[~, IA, ~] = uniquetol([IC1, IC2], 1e-5, 'ByRows', true, 'OutputAllIndices', true);

n = size(A,2);

% Count how many rows belong to each unique cluster
counts = cellfun(@numel, IA);
numClusters = length(IA);
C = zeros(numClusters,2*n);
for j = 1:numClusters
    C(j,:) = M(IA{j}(1), :);
end

flaggedIndices = [];
visited = false(numClusters, 1);

for i = 1:numClusters
    if visited(i), continue; end
    
    currentCenter = M(IA{i}(1), :);
    currentCount = counts(i);
    
    % Check if cluster is Real or Complex
    isReal = all(abs(currentCenter(n+1:2*n)) < imagtol);
    
    if isReal
        % Logic: Find real repeated rows
        if currentCount > 1
            flaggedIndices = [flaggedIndices; IA{i}];
        end
        visited(i) = true;
    else
        % Logic: Find complex clusters without matching conjugate count
        conjCenter = [currentCenter(1:n) -currentCenter(n+1:2*n)];
        
        % Find the cluster index that matches the conjugate center
        isMatch = false;
        for j = 1:numClusters
            if norm(conjCenter - C(j,:))<tol*(1+norm(conjCenter))
                conjIdx = j;
                isMatch = true;
            end
        end
        
        if isMatch
            conjCount = counts(conjIdx);
            if currentCount ~= conjCount
                % Flag the one that has more rows
                targetIdx = i;
                if conjCount > currentCount
                    targetIdx = conjIdx;
                end
                flaggedIndices = [flaggedIndices; IA{targetIdx}];
            end
            visited(i) = true;
            visited(conjIdx) = true;
        else
            % No conjugate cluster found at all
            flaggedIndices = [flaggedIndices; IA{i}];
            visited(i) = true;
        end
    end
end

ind = unique(flaggedIndices);