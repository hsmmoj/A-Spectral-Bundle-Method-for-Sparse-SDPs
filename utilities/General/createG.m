function G = createG(P, clique, dim)
    G = {};
    V = {};
    for i=1:length(P)
        V{i} = LHSExp(P{i}, clique{i}, dim);
    end
    
    for i=1:length(V)
        for j=1:length(V)
            dummyG = V{i}' * V{j};
            r = length(dummyG);
            dummyG2 = zeros(r^2);
            for n=1:r
                for m=1:r
                    dummyG2((m-1)*r+1:m*r,(n-1)*r+1:n*r) = dummyG(:,n) * dummyG(m,:);
                end
            end
             G{(i-1)*r + j} = dummyG2;
        end
    end
end