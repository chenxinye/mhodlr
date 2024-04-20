function T = htriu(A)
    T = A;
    if isempty(A.D)
        T.D = triu(A.D);
    else
        T.U2 = zeros(size(A.U2, 1), 0);
        T.V1 = zeros(size(A.V1, 1), 0);
        T.A11 = htriu(A.A11);
        T.A22 = htriu(A.A22);
    end
end