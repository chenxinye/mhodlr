function [L, U] = hlu(H, threshold)
    if ~issquare(H)
        error('LU factorization is only applied to a square HODLR matrix.');
    end
    
    [m, n, m1, m2, n1, n2] = hsize(H);
    
    if isempty(H.D)
        [L11, U11] = hlu(H.A11, threshold);
        % U12 = mldivide(L11, H.U1 * H.V2); % L11 U12 = A12 = H.U1 * H.V2
        % L21 = mrdivide(H.U2 * H.V1, U11); % L21 U11 = A21 = H.U2 * H.V1
        U12 = mldivide(L11, H.U1); 
        L21 = mrdivide(H.V1, U11); 

        X = -H.U2 * (L21 * U12); 
        NH = hrank_update(H.A22, X, H.V2, threshold);
        [L22, U22] = hlu(NH, threshold); % lu(hadd(H.A22, L21 * U12, '-'));
        
        U12 = U12 * H.V2;
        L21 = H.U2 * L21;

        L = blkdiag(L11, L22);
        L(size(L11, 1)+1:end, 1:size(L11, 2)) = L21;
        L = sparse(L);
        U = blkdiag(U11, U22);
        U(1:size(U11, 1), size(U11, 2)+1:end) = U12;
        U = sparse(U);
    else
        [L, U] = lu(H.D);
    end
end 