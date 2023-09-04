function bigMatrix = matrixExp(matrix, clique, size)
bigMatrix = zeros(size,size);
m = length(clique);
linearIndices = sub2ind([size size],repmat(clique,[1,m]), reshape(repmat(clique,[m,1]), [1,m*m]));
bigMatrix(linearIndices) = matrix;

end