% function L =createL(P, W, clique, dim,Paras)
% %     L = {};
% %     p = length(clique);
% %     for i=1:p
% %         for j=1:p
% %             M = 
% %             L{(i-1)*p + j} = 2 * M';
% %         end
% %     end
% 
%     L = {};
%     %V = {};
%     overlaps = Paras.overlaps;
%     z = 1;
%     ind = {};
%     Bigind = {};
%     for i=1:length(clique)
%         for j=1:length(clique)
%             [~,H1] = size(clique{i});
%             [~,H2] = size(clique{j});
%             pairs = overlaps{z};
%             for k1 = 1:height(pairs)
%                  ind{end+1} = [(pairs(k1,1)-1) * H1 + pairs(:,1), (pairs(k1,2)-1) * H2 + pairs(:,2)];
%             end
%             Bigind{z} = vertcat(ind{:});
%             ind = {};
%             z = z + 1;
%         end 
%     end
% 
%     kronP = {};
%     for i =1:length(clique)
%         kronP{i} = kron(P{i},P{i});
%     end
% 
%     z = 1;
%     for i=1:length(clique)
%         for j=1:length(clique)
%             if i == j 
%                P1 = kronP{j};
%                P2 = vec(W{i});
%                L{z} = P2'*P1*2;
%             elseif Bigind{z} == 0
%                L{z} = zeros(1,Paras.MaxCols^2);
%             else
%                P1 = kronP{j};
%                %P2 = kronP{j};
%                P2 = vec(W{i});
%                L{z} = zeros(1,Paras.MaxCols^2);
%                indpair = Bigind{z};
%                for k = 1:length(indpair) 
%                    L{z} = L{z} + P1(indpair(k,1),:)*P2(indpair(k,2))*2;
%                end
%             end
%             z = z + 1;
%         end
%     end
% end
% 














function L =createL(P, W, clique, dim,Paras)
    L = {};
    p = length(clique);
    for i=1:p
        for j=1:p
            M = vec(P{j}'*subMatrixExt(matrixExp(W{i},clique{i},dim), clique{j})*P{j});
            L{(i-1)*p + j} = 2 * M';
        end
    end
end