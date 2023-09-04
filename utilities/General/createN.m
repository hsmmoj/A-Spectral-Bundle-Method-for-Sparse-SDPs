function N = createN(A)
     m = size(A,1);
     n = size(A,3);
    dummy = reshape(A, [m^2,n])';
    N = dummy * dummy';
end