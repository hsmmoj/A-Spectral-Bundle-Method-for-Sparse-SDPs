function [X_star, count, Out] = BundleDecomposedSDP(A, B, C, opts)
    % Preprocess the problem data
    [cliques, W, omega, P, numOfConstraints, numOfCliques, overlaps, subA, subC, N, invQ33] = preproccesing(A, B, C, opts);
    count = 0;
    alpha = [];
    Obj = [];
    DualAffine = [];
    Gap = [];
    
    % Define parameters for the optimization problem
    Paras.MaxCols = opts.r_c + opts.r_p;
    Paras.NumOfVar = Paras.MaxCols^2 + 1;
    Paras.numOfCliques = numOfCliques;
    Paras.numOfConstraints = numOfConstraints;
    Paras.overlaps = overlaps;
    Paras.subA = subA;
    Paras.subC = subC;
    Paras.N = N;
    Paras.invQ33 = invQ33;

    % Calculate indices for symmetric variables
    for i = 1:numOfCliques
        [xIndSym, ~, xIndOffDiag, ~, ~, xIndOffDiagCounter] = SymmetricIndices(Paras.MaxCols, false);
        Paras.IndicesPSD{i} = xIndSym;
        Paras.IndOffDiagPSD{i} = xIndOffDiag;
        Paras.IndOffDiagCounter{i} = xIndOffDiagCounter;
    end

    [XIndSym, ~, XIndOffDiag, ~, ~, XIndOffDiagCounter] = SymmetricIndices(opts.dim, false);
    Paras.XIndSym = XIndSym;
    Paras.XIndOffDiag = XIndOffDiag;
    Paras.XIndOffDiagCounter = XIndOffDiagCounter;

    for iter = 1:opts.maxIter
        % Solve the Master Problem
        [X_star, gamma_val, W_val, S_val, y_val, dualAffine, gap] = MasterProblem(omega, A, B, C, W, P, cliques, numOfCliques, opts, Paras);
        
        % Check the termination criteria
        [decent, stop, realValAtOmega, optsNew] = criteria(omega, X_star, C, numOfCliques, cliques, W_val, opts, iter);
        opts = optsNew;
        
        % Store optimization results
        Obj(iter) = realValAtOmega;
        alpha(iter) = opts.alpha;
        DualAffine(iter) = dualAffine;
        Gap(iter) = gap;
        
        % Update the current solution if needed
        if iter == 1
            omega = X_star;
        else
            if stop && iter > 1
                break;
            end
            if decent
                omega = X_star;
                count = count + 1;
            end
        end
        
        % Compute next parameters W and P
        [W, P] = computeNextPar(X_star, W, P, S_val, gamma_val, numOfCliques, cliques, W_val, opts);
        
        % Display iteration progress
    fprintf('Iteration: %5d\n', iter);
    fprintf('Objective Value: %7.5e\n', realValAtOmega);
    fprintf('Alpha Value: %7.2e\n', opts.alpha);
    fprintf('Dual Affine: %7.5e\n', dualAffine);
    fprintf('Dual Gap: %7.5e\n', gap);
    fprintf('-------------------------\n');  % Separator line
    end

    % Store the optimization results in the 'Out' structure
    Out.Obj = Obj;
    Out.alpha = alpha;
    Out.omega = omega;
    Out.dualAffine = DualAffine;
    Out.gap = Gap;
end
