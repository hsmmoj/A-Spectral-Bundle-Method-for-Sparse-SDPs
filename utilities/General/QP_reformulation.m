function [Q, q] = QP_reformulation(A, B, C, W, P, r_p, r_c, cliques, omega, alpha, Paras, opts)
    % Compute the dimensions and initialize Q matrix
    r = r_c + r_p;
    Numcliques = length(cliques);
    NumConstraints = size(A, 3);
    dimQ = Numcliques + (r_p + r_c)^2 * Numcliques + NumConstraints;
    dim = size(A, 1);
    Q = zeros(dimQ);

    % Fill the Q matrix with components
    indH = sub2ind([dimQ dimQ], repmat(linspace(1, Numcliques, Numcliques), [1, Numcliques]), reshape(repmat(linspace(1, Numcliques, Numcliques), [Numcliques, 1]), [1, Numcliques^2]));
    Q(indH) = cell2mat(createH(W, cliques, dim, Paras));

    indL = sub2ind([dimQ dimQ], repmat(linspace(1, Numcliques, Numcliques), [1, (r^2) * Numcliques]), reshape(repmat(linspace(Numcliques + 1, r^2 * Numcliques + Numcliques, (r^2) * Numcliques), [Numcliques, 1]), [1, (r^2) * Numcliques^2]));
    L = mat(cell2mat(createL(P, W, cliques, dim, Paras)));
    Q(indL) = reshape(L, [r^2 * Numcliques, Numcliques])';

    G = createG2(P, cliques, dim, Paras);
    Grow = 0;
    Gcol = 0;
    for i = 1:length(G)
        indG = sub2ind([dimQ dimQ], repmat(linspace(Numcliques + 1 + Grow * r^2, r^2 + Numcliques + Grow * r^2, r^2), [1, (r^2)]), reshape(repmat(linspace(Numcliques + 1 + r^2 * Gcol, r^2 + Numcliques + r^2 * Gcol, (r^2)), [r^2, 1]), [1, (r^4)]));
        Q(indG) = G{i};
        Gcol = Gcol + 1;
        if mod(i, Numcliques) == 0
            Grow = Grow + 1;
            Gcol = 0;
        end
    end

    indF = sub2ind([dimQ dimQ], repmat(linspace(1, Numcliques, Numcliques), ...
        [1, NumConstraints]), reshape(repmat(linspace(Numcliques + r^2 * ...
        Numcliques + 1, NumConstraints + Numcliques + r^2 * Numcliques, ...
        NumConstraints), [Numcliques, 1]), [1, NumConstraints * Numcliques]));
    Q(indF) = reshape(cell2mat(createF(A, W, cliques, dim, Paras)), ...
        [NumConstraints, Numcliques])';

    indR = sub2ind([dimQ dimQ], repmat(linspace(Numcliques + 1, Numcliques + r^2 * Numcliques, r^2 * Numcliques), [1, NumConstraints]), reshape(repmat(linspace(Numcliques + r^2 * Numcliques + 1, NumConstraints + Numcliques + r^2 * Numcliques, NumConstraints), [r^2 * Numcliques, 1]), [1, NumConstraints * Numcliques * r^2]));
    RDummy = createR(A, P, cliques, Paras);
    R = [];
    for i = 1:length(RDummy)
        R = [R; RDummy{i}];
    end
    Q(indR) = R;

    indN = sub2ind([dimQ dimQ], repmat(linspace(Numcliques + Numcliques * r^2 + 1, NumConstraints + Numcliques + Numcliques * r^2, NumConstraints), [1, NumConstraints]), reshape(repmat(linspace(Numcliques + r^2 * Numcliques + 1, NumConstraints + Numcliques + r^2 * Numcliques, NumConstraints), [NumConstraints, 1]), [1, NumConstraints^2]));
    Q(indN) = Paras.N;

    % Ensure Q is symmetric
    Q = 0.5 * (Q + Q');

    % Initialize q vector
    q = zeros(dimQ, 1);

    % Fill the q vector
    J = createJ(W, C, omega, cliques, alpha);
    q(1:Numcliques) = cell2mat(J);
    q(Numcliques + 1:r^2 * Numcliques + Numcliques) = createM(P, C, omega, cliques, alpha, Paras);
    q(Numcliques + r^2 * Numcliques + 1:end) = createP(A, B, C, omega, alpha, opts);
end
