function T = htriu(H)
    T = H;
    if isempty(H.D)
        T.D = triu(H.D);
    else
        T.U2 = zeros(size(H.U2, 1), 0);
        T.V1 = zeros(size(H.V1, 1), 0);
        T.A11 = htriu(H.A11);
        T.A22 = htriu(H.A22);
    end
end