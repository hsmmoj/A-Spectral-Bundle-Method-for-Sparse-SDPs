function M = createM(P, C, omega, cliques, alpha, Paras)
    M = [];
    %M = zeros(1,length(cliques)*Paras.MaxCols^2);
    for i=1:length(cliques)
        %M((i-1)*Paras.MaxCols^2+1:i*Paras.MaxCols^2)=kron(P{i}',P{i}')*(Paras.subC{i}*(-2) + 2*alpha*subMatrixExt_Vec(omega,cliques{i}));
        M = [M vec(P{i}'*subMatrixExt(-2*C + 2*alpha*omega,cliques{i})*P{i})'];
    end
end