function H = createH(W, clique, dim,Paras)
% 
    H = {};
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

%     kronP = {};
%     for i =1:length(clique)
%         kronP{i} = kron(P{i},P{i});
%     end

    z = 1;
    for i=1:length(clique)
        for j=1:length(clique)
             P1 = vec(W{i});
             P2 = vec(W{j});
            if i == j 
               H{z} = P1'*P2;%eye(Paras.MaxCols^2);
            elseif Bigind{z} == 0
               H{z} = 0;%zeros(Paras.MaxCols^2);
            else
%                P1 = kronP{i};
%                P2 = kronP{j};
              
               %H{z} = zeros(Paras.MaxCols^2);
               H{z} = 0;
               indpair = Bigind{z};
               for k = 1:length(indpair) 
                   %H{z} = H{z} + P1(indpair(k,1),:)'*P2(indpair(k,2),:);
                   H{z} = H{z} + P1(indpair(k,1))*P2(indpair(k,2));
               end
            end
            z = z + 1;
        end
    end










%     H = {};
%     r = length(W);
%     for i=1:r
%         for j=1:r
%             dummy1 = matrixExp(W{i}, clique{i}, dim);
%             dummy2 = matrixExp(W{j}, clique{j}, dim);
%             H{(i-1)*r + j} = trace(dummy1 * dummy2');
%         end
%     end
end          