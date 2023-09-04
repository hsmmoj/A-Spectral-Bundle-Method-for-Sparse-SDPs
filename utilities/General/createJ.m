function J = createJ(W, C, omega, cliques, alpha)
    J = {};
    for i=1:length(cliques)
        J{i} = trace(W{i}' * (-2 * subMatrixExt(C, cliques{i}) + 2 * alpha * subMatrixExt(omega, cliques{i})));
    end
end