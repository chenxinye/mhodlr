function T = htril(H)
%{
    Return the lower triangular part of HODLR matrix in HODLR format.
%}
    T = H;
    if isempty(H.D)
        T.D = tril(H.D);
    else
        T.U1 = zeros(size(H.U1, 1), 0);
        T.V2 = zeros(size(H.V2, 1), 0);
        T.A11 = htril(H.A11);
        T.A22 = htril(H.A22);
    end
end