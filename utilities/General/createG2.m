function G = createG2(P, clique, dim, Paras)
    G = {};
    %V = {};
    overlaps = Paras.overlaps;
    z = 1;
    ind = {};
    Bigind = {};
    for i=1:length(clique)
        for j=1:length(clique)
            [~,H1] = size(clique{i});
            [~,H2] = size(clique{j});
            pairs = overlaps{z};
            for k1 = 1:height(pairs)
                 ind{end+1} = [(pairs(k1,1)-1) * H1 + pairs(:,1), (pairs(k1,2)-1) * H2 + pairs(:,2)];
            end
            Bigind{z} = vertcat(ind{:});
            ind = {};
            z = z + 1;
        end 
    end

    kronP = {};
    for i =1:length(clique)
        kronP{i} = kron(P{i},P{i});
    end

    z = 1;
    for i=1:length(clique)
        for j=1:length(clique)
            if i == j 
               G{z} = eye(Paras.MaxCols^2);
            elseif Bigind{z} == 0
               G{z} = zeros(Paras.MaxCols^2);
            else
               P1 = kronP{i};
               P2 = kronP{j};
               G{z} = zeros(Paras.MaxCols^2);
               indpair = Bigind{z};
               for k = 1:length(indpair) 
                   G{z} = G{z} + P1(indpair(k,1),:)'*P2(indpair(k,2),:);
               end
            end
            z = z + 1;
        end
    end



%     for i=1:length(P)
%         V{i} = LHSExp(P{i}, clique{i}, dim);
%     end
%     
%     for i=1:length(V)
%         for j=1:length(V)
%              G{(i-1)*length(V) + j} = kron(V{i}',V{i}') * kron(V{j},V{j}) ;
%         end
%     end
end