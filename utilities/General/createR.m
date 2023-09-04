function R = createR(A, P, clique, Paras)
    R = {};
    for i=1:length(clique)
        R{i} = 2 * kron(P{i}', P{i}') * Paras.subA{i};
    end  
end 



% function R = createR(A, P, clique,Paras)
%     R = {};
%     for i=1:length(clique)
%         dummyR = [];
%         for j=1:size(A,3)
%             dummyR = [dummyR kron(P{i}', P{i}')*vec(subMatrixExt(A(:,:,j), clique{i}))];
%         end
%         R{i} = 2 * dummyR;
%     end  
% end 

% 
% function R = createR(A, P, clique, Paras)
%     R = {};
%     At = vec(A);
%     for i=1:length(clique)
%         dummyR = zeros(width(P{i})^2 ,size(A,3));
%         for j=1:size(A,3)
%             dummyR(:,j) = kron(P{i}', P{i}')*vec(subMatrixExt(A(:,:,j), clique{i}));
%         end
%         R{i} = 2 * dummyR;
%     end  
% end 
