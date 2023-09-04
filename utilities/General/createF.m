function F = createF(A, W, clique, dim, Paras)
    F = {};
    for i=1:length(clique)
        F{i} = 2 *  vec(W{i})'* Paras.subA{i};
    end

end