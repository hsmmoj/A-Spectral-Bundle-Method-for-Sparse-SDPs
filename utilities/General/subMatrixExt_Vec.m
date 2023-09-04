function subMatrix = subMatrixExt_Vec(matrix, clique)
m1 = length(clique);
linearIndices = sub2ind(size(matrix),repmat(clique,[1,m1]), reshape(repmat(clique,[m1,1]), [1,m1*m1]));
%subMatrix = zeros(length(clique));
subMatrix=matrix(linearIndices)';
%subMatrix = mat(matrix(linearIndices));
end