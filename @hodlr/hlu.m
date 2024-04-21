function [L, U] = hlu(A)
    if ~issquare(A)
        error('LU factorization is only applied to a square HODLR matrix.');
    end
    
    [m, n, m1, m2, n1, n2] = hsize(A);
    
    if isempty(A.D)
        [L11, U11] = hlu(A.A11);
        U12 = mldivide(L11, A.U1 * A.V2); % L11 U12 = A12
        L21 = mrdivide(A.U2 * A.V1, U11); % L21 U11 = A21
        [L22, U22] = lu(hadd(A.A22, L21 * U12, '-'));
        
        L = blkdiag(L11, L22);
        L(size(L11, 1)+1:end, 1:size(L11, 2)) = L21;
        L = sparse(L);
        U = blkdiag(U11, U22);
        U(1:size(U11, 1), size(U11, 2)+1:end) = U12;
        U = sparse(U);
    else
        [L, U] = lu(A.D);
    end
end 