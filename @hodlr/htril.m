function T = htril(A)
    T = A;
    if isempty(A.D)
        T.D = tril(A.D);
    else
        T.U1 = zeros(size(A.U1, 1), 0);
        T.V2 = zeros(size(A.V2, 1), 0);
        T.A11 = htril(A.A11);
        T.A22 = htril(A.A22);
    end
end