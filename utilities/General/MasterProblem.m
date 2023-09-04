function [X_star, gamma_val, W_val, S_val, y_val, DualFeasibility, gap] = MasterProblem(omega, A, B, C, W, P, cliques, numOfCliques, opts, Paras)
    % Preprocess problem data
    b_sdp = B;
    c_sdp = vec(C);

    % Reformulate the problem into a Quadratic Programming
    [Q, q] = QP_reformulation(A, B, C, W, P, opts.r_p, opts.r_c, cliques, omega, opts.alpha, Paras, opts);
    
    % Extract matrices and vectors from Q and q
    r = opts.r_c + opts.r_p;
    Q11 = Q(1:numOfCliques, 1:numOfCliques);
    Q12 = Q(1:numOfCliques, numOfCliques + 1:numOfCliques + r^2 * numOfCliques);
    Q13 = Q(1:numOfCliques, numOfCliques + r^2 * numOfCliques + 1:end);
    Q22 = Q(numOfCliques + 1:numOfCliques + r^2 * numOfCliques, numOfCliques + 1:numOfCliques + r^2 * numOfCliques);
    Q23 = Q(numOfCliques + 1:numOfCliques + r^2 * numOfCliques, numOfCliques + r^2 * numOfCliques + 1:end);
    Q33 = Q(numOfCliques + r^2 * numOfCliques + 1:end, numOfCliques + r^2 * numOfCliques + 1:end);
    q1 = q(1:numOfCliques);
    q2 = q(numOfCliques + 1:numOfCliques + r^2 * numOfCliques);
    q3 = q(numOfCliques + 1 + r^2 * numOfCliques:end);

    % Compute matrices M
    M11 = Q11 - Q13 * Paras.invQ33 * Q13';
    M22 = Q22 - Q23 * Paras.invQ33 * Q23';
    M12 = Q12 - Q13 * Paras.invQ33 * Q23';
    m1 = q1 - Q13 * Paras.invQ33 * q3;
    m2 = q2 - Q23 * Paras.invQ33 * q3;
    M = [M11, M12; M12', M22];
    M = (M + M') / 2;

    [eig_vec, eig_val] = eig(M);
    % Ensure eigenvalues are non-negative (numerical stability)
    eig_val(eig_val < 0) = 0;
    M05 = eig_vec * sqrt(eig_val) * eig_vec.';
    M05 = (M05 + M05.') / 2;

    % Ensure symmetry for specific elements
    start = 0;
    for i = 1:numOfCliques
        M05(:, numOfCliques + start + [Paras.IndOffDiagPSD{i}; Paras.IndOffDiagCounter{i}]) = ...
            1/2 * (M05(:, numOfCliques + start + [Paras.IndOffDiagPSD{i}; Paras.IndOffDiagCounter{i}]) + ...
            M05(:, numOfCliques + start + [Paras.IndOffDiagCounter{i}; Paras.IndOffDiagPSD{i}]));
        start = start + Paras.MaxCols^2;
    end

    % Prepare constraint matrices for the Sedumi solver
    [height_M05, ~] = size(M05);
    At = zeros(height_M05 + numOfCliques + 1, (2) * numOfCliques + 2 + height_M05 + numOfCliques * Paras.MaxCols^2);
    At(1:height_M05, numOfCliques + 1:2 * numOfCliques) = M05(:, 1:numOfCliques);
    At(1:height_M05, numOfCliques * (1 + 1) + (2 + height_M05) + 1:end) = M05(:, numOfCliques + 1:end);
    At(1:height_M05, numOfCliques * (1 + 1) + 2 + 1:numOfCliques * (1 + 1) + 2 + height_M05) = -eye(height_M05);
    for i = 1:numOfCliques
        At(height_M05 + (i - 1) + 1, i) = 1; % a
        At(height_M05 + (i - 1) + 1, numOfCliques + i) = 1; % gamma
        At(height_M05 + (i - 1) + 1, 2 * numOfCliques + 2 + height_M05 + (i - 1) * Paras.MaxCols^2 + 1:2 * numOfCliques + 2 + height_M05 + (i) * Paras.MaxCols^2) = reshape(eye(Paras.MaxCols), 1, []); % s
    end
    At(end, 2 * numOfCliques + 1 + 1) = 1;

    % Prepare right-hand side vector b and the cost vector c
    b = zeros(height_M05 + numOfCliques + 1, 1);
    b(height_M05 + 1:end - 1) = ones(numOfCliques, 1) * opts.rho;
    b(end) = 0.5;
    c = [zeros(numOfCliques, 1); m1; 1; 0; zeros(height_M05, 1); m2];

    % Define optimization problem dimensions
    K.s = ones(1, numOfCliques) * Paras.MaxCols;
    K.r = 2 + height_M05;
    K.l = 2 * numOfCliques;

    % Call MOSEK to solve the optimization problem
    prob1 = SedumiToMosek_Latest(At, b, c, K);
    [~, res1] = mosekopt('minimize echo(0)', prob1);
    status = res1.sol.itr.prosta;
    
    % Check if the solution is feasible
    if ~strcmp(status, 'PRIMAL_AND_DUAL_FEASIBLE')
        warning('Infeasible solution');
        X_star = [];
        gamma_val = [];
        W_val = {};
        S_val = {};
        y_val = [];
        DualFeasibility = NaN;
        gap = NaN;
        return;
    end

    % Extract solution variables
    xlinear = res1.sol.itr.xx;
    xPSD = res1.sol.itr.barx;

    Gammastar = xlinear(numOfCliques + 1:2 * numOfCliques);
    VecSstar = [];
    for i = 1:numOfCliques
        Sstar{i} = zeros(Paras.MaxCols);
        Sstar{i}(Paras.IndicesPSD{i}) = xPSD((i - 1) * (Paras.MaxCols + 1) * Paras.MaxCols / 2 + 1:(i) * (Paras.MaxCols + 1) * Paras.MaxCols / 2);
        Sstar{i}(Paras.IndOffDiagCounter{i}) = Sstar{i}(Paras.IndOffDiagPSD{i});
        VecSstar = [VecSstar; vec(Sstar{i})];
        Wstar{i} = Gammastar(i) * W{i} + P{i} * Sstar{i} * P{i}.';
    end

    y = Paras.invQ33 * (-q3 / 2 - Q13' * Gammastar - Q23' * VecSstar);

    % Improve numerical stability
    DualAffine = 0;
    for i = 1:numOfCliques
        DualAffine = DualAffine + matrixExp(Wstar{i}, cliques{i}, opts.dim);
    end
    DualAffine = vec(DualAffine);
    DualAffine = DualAffine - c_sdp + opts.At * y;
    DualAffine([Paras.XIndOffDiag, Paras.XIndOffDiagCounter]) = ...
        1/2 * (DualAffine([Paras.XIndOffDiag, Paras.XIndOffDiagCounter]) + ...
        DualAffine([Paras.XIndOffDiagCounter, Paras.XIndOffDiag]));

    % Compute the next solution and evaluate dual feasibility and gap
    X_next = vec(omega) + (DualAffine) / opts.alpha;
    DualFeasibility = norm(DualAffine, 'fro')^2 / (1 + norm(C, 'fro'));
    gap = abs(b_sdp' * y - c_sdp' * X_next) / (1 + abs(b_sdp' * y) + abs(c_sdp' * X_next));
    X_star = mat(X_next);
    gamma_val = Gammastar;
    for i = 1:numOfCliques
        W_val{i} = Wstar{i};
        S_val{i} = Sstar{i};
    end
    y_val = y;
end
