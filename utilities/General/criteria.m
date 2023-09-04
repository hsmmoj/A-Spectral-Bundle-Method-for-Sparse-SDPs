function [decent, stop, realValAtOmega, opts] = criteria(omega, X_star, C, numOfCliques, cliques, W_val, opts, iter)
    % Compute the real values at omega and X_star
    realValAtOmega = trace(C' * omega);
    realValAtXstar = trace(C' * X_star);
    appValAtXstar = realValAtXstar;

    % Compute real values for each clique and update the overall real value
    for i = 1:numOfCliques
        smallX = subMatrixExt(X_star, cliques{i});
        [eigVecX, eigValX] = eig(-smallX, 'vector');
        smallOmega = subMatrixExt(omega, cliques{i});
        [eigVecOmega, eigValOmega] = eig(-smallOmega, 'vector');
        realValAtOmega = realValAtOmega + opts.rho * max(0, max(eigValOmega));
        realValAtXstar = realValAtXstar + opts.rho * max(0, max(eigValX));
        appValAtXstar = appValAtXstar - trace(W_val{i}' * smallX);
    end

    stop = 0;
    decent = 0;
    null = 0;

    % Check the criteria for deciding if the solution is decent
    if opts.beta * (realValAtOmega - appValAtXstar) <= realValAtOmega - realValAtXstar
        decent = 1;
        
        % Check if the dynamic alpha option is enabled
        if (opts.mu * (realValAtOmega - appValAtXstar) <= realValAtOmega - realValAtXstar) && opts.dynamicAlpha
            opts.alpha = max(10^-5, opts.alpha /1.1);
        end

        % Additional alpha modification conditions
%         if iter > 100
%             opts.alpha = 0.0001;
%         elseif null == 10
%             opts.alpha = 0.1;
%             null = 0;
%         end

        % % Check if the dynamic alpha option is enabled
        % if opts.ml * (realValAtOmega - appValAtXstar) >= (realValAtOmega - realValAtXstar) && opts.dynamicAlpha
        %     opts.alpha = min(opts.alpha * 1.01, 1);
        % end
        % null = null + 1;
    end

    % Check if the stopping criterion is met
    if (realValAtOmega - appValAtXstar) <= opts.eps
        stop = 1;
    end
end