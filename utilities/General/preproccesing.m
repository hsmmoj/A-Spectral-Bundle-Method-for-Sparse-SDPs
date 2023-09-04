function [cliques, W, omega, P, numOfConstraints, numOfCliques, overlaps, subA, subC, N, invQ33] = preprocessing(A, B, C, opts)
    import cdcs_utils.*;
    
    % Extract necessary data
    numOfConstraints = size(B, 1);
    K.l = 0;
    K.q = 0;
    K.f = 0;
    K.s = sqrt(size(opts.At, 1));
    [~, ~, Ats, Cs] = svecData(opts.At, opts.c, K);
    temp = chordalDecomposition(Ats, Cs, K);
    cliques = temp{1}.Set;
    
    % Adjust clique indices
    for i = 1:length(cliques)
        cliques{i} = cliques{i} - 1;
    end
    
    % Initialize variables
    overlaps = {};
    subA = {};
    subC = {};
    
    % Compute clique overlaps
    for i = 1:length(cliques)
        for j = 1:length(cliques)
            [~, ia, ib] = intersect(cliques{i}, cliques{j});
            overlaps{end + 1} = [ia, ib];
        end
    end

    numOfCliques = length(cliques);
    Ec = {};
    
    % Construct Ec
    for k = 1:length(cliques)
        dummy = zeros(length(cliques{k}), opts.dim);
        for i = 1:length(cliques{k})
            dummy(i, cliques{k}(i) + 1) = 1;
        end
        Ec = [Ec; dummy];
    end
    
    % Initialize Omega, P, and W
    omega = eye(opts.dim);
    P = {};
    W = {};
    V = {};
    
    % Compute P, W, and submatrices for each clique
    for i = 1:length(cliques)
        cliques{i} = cliques{i} + 1;
        omegaSmall = subMatrixExt(omega, cliques{i});
        [eigVec, eigVal] = eig(-omegaSmall);
        [d, ind] = sort(diag(eigVal), 'descend');
        Vs = eigVec(:, ind);
        V{i} = Vs(:, 1:opts.r_c);
        P{i} = Vs(:, 1:(opts.r_p + opts.r_c));
        W{i} = Vs(:, 1) * Vs(:, 1)';
        
        % Compute submatrices for A and C
        tempA = zeros(length(cliques{i})^2, size(A, 3));
        for j = 1:size(A, 3)
            tempA(:, j) = vec(subMatrixExt(A(:, :, j), cliques{i}));
        end
        subA{i} = tempA;
        subC{i} = vec(subMatrixExt(C, cliques{i}));
    end
    
    % Compute N and invQ33
    N = createN(A);
    U = chol(Ats.' * Ats);
    inv_U = inv(U);
    invQ33 = inv_U * inv_U.';
end
