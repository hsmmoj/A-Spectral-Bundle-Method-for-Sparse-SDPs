function [newW, newP, newV] = computeNextPar(X_star,oldW, oldP, S_val, gamma_val,numOfCliques,cliques,W_val,opts)
    newW = {};
    newP = {};
    newV = {};
    for i=1:numOfCliques
        smallX = subMatrixExt(X_star, cliques{i});
        [eigVec, eigVal]    = eig(-smallX);
        [d,ind]             = sort(diag(eigVal), 'descend');
        sortedEigVecX       = eigVec(:,ind);
        newV{i} = sortedEigVecX(:,1:opts.r_c);
        if opts.r_p ~= 0
            [eigVec, eigVal]    = eig(S_val{i});
            [sortedEigValS,ind] = sort(diag(eigVal), 'descend');
            sortedEigVecS       = orth(eigVec(:,ind));
            Q1 = sortedEigVecS(:,1:opts.r_p);
            Q2 = sortedEigVecS(:,opts.r_p+1:end);
            Sigma2 = diag(sortedEigValS(opts.r_p+1:end));
            newP{i}             = orth([newV{i} oldP{i}*Q1]);
    
            if gamma_val > 0
                newW{i}         =  ( gamma_val(i) * oldW{i} + newP{i}*Q2*Sigma2*Q2'*newP{i}') /(gamma_val(i) + trace(Sigma2));
            else
               newW{i}         = oldW{i};
           end
        else
            newP{i} = newV{i};
            newW{i} = W_val{i} / trace(W_val{i});
        end
    end
end