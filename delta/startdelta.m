function delta = startdelta(A,version)
    %STARTDELTA   Delta matrices for start problem.
    %   delta = STARTDELTA(A,m) provides the delta matrices for an 
    %   m-parameter start problem with constant coefficient matrix A.
    
    % version: 
    %   1 - we build perms(1:m) and reduce from all possible permutations
    %   2 - we build possible permutations from nonzero indices 

    if nargin<2
        version=2;
    end

    % Retrieve the matrix sizes:
    [k, l] = size(A);
    m = k - l + 1;

    % Determine the left and right compression matrices:
    L = left_compression_mat(k, m);
    T = right_compression_matrix(l, m);

    delta = cell(m + 1,1);
    % we compute all permutations only once
    if version==1
        allperms = perms(1:m);
    else
        allperms = [];
    end
    for i = 1:m + 1
        % Compute all the nonzero permutations:
        plist = [];
        for j= 1:size(L, 1)
            % Compute all the nonzero permutations for this row:
            permutation = permutations(L(j,:), m, l, i - 1, version, allperms);

            for idp = 1:size(permutation, 1)
                % Determine the elements in A that are selected:
                if i == 1
                    a = 1;
                    ida = 0;
                else
                    idrow = L(j, :);
                    a = A(idrow(permutation(idp, :) == (i - 1)), :);
                    ida = find(a).';
                    a = a(ida).';
                end

                % Build list with idcol, permutations, coefficients, and 
                % positions of these coefficients:
                p = [j*ones(length(a), 1) kron(permutation(idp, :),...
                     ones(length(a), 1)) a ida];
                plist = [plist; p];
            end
        end

        % Pre-allocate the rows, columns, and elements:
        % (rows and elements up to sign are already in plist.)
        row = plist(:, 1);
        col = zeros(size(plist, 1), 1);
        ell = plist(:, end - 1);
        for j = 1:size(plist, 1)
            % Retrieve row multi-index and permutation:
            idrow = L(plist(j, 1), :);
            ida = plist(j, end);
            permutation = plist(j, 2:end - 2);
            
            % Compute Levi-Civita symbol:
            P = speye(m);
            sign = det(P(:, permutation));
            
            % Find row index, column index, and element:
            idcol = subkron(permutation, idrow, i - 1, ida, l, m);
            col(j) = idcol;
            ell(j) = sign*ell(j);
        end

        % Construct delta matrix:
        delta{i} = sparse(row, col, ell, size(L,1), l^m)*T;
    end
end 

%% Define auxiliary functions:
function [idcol] = subkron(permutation, idrow, iddelta, ida, l, m)
    %SUBKRON   Perform Kronecker operation.
    %   [idcol] = SUBKRON(permutation, idrow, iddelta, ida, l, m) finds the
    %   column index of the Kronecker operation for a specific permutation,
    %   given the row multi-index. It works for the pre-specified delta
    %   matrix (with iddelta) and special identity coefficient matrices.
    %
    %   The position of the coefficient is given by ida.
    
    if iddelta == 0
        idcol = idrow - permutation + 1;
    else
        % Compute idcol (expanded) for each selected element of A:
        idcol = idrow - permutation + 1;
        idcol(permutation == iddelta) = ida;
    end

    % Change idcol (expanded) into a single column index:
    idcol = sum((idcol - 1).*l.^(m - 1:-1:0), 2) + 1;
end

function p = permutations(idrow, m, l, iddelta, version, allperms)
    %PERMUTATIONS   Nonzero-permutation.
    %   p = PERMUTATIONS(idx, m, l, iddelta) returns all
    %   nonzero-permutations in the Kronecker operation for a given row
    %   multi-index.

    % Determine all permutations:
    p = allperms;

    % Remove zero-permutations for each index (version 1):
    G = cell(1,m);
    for i = 1:m
        v = idrow(i) >= 1:m;
        w = idrow(i) <= l:l+m-1;
        z = all([v; w]);
        if iddelta ~= 0
            z(iddelta) = 1;
        end
        numlist = 1:m;
        G{i} = numlist(z);
        if version==1
            mask = ismember(p(:,i), numlist(z));
            p = p(mask,:);
        end
    end
    if version == 2
        % Build pertmutations from lists of possible indices
        p = fast_perms(G);
    end
end

function p = fast_perms(z)
    
    k = length(z);
    p = z{1}';
    for j = 2:k
        newp = [];
        for q = 1:length(z{j})
            ismb = ismember(p, z{j}(q));
            if size(ismb,2)>1
                ismb = any(ismb');
            end
            mask = not(ismb);
            if any(mask)
                newp = [newp; p(mask,:) z{j}(q)*ones(sum(mask),1)];
            end
        end
        p = newp;
    end
end