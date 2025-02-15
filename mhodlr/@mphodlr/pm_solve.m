function x = pm_solve(H, b)
%{
    A recursive way to solve linear system Ax = b with HODLR matrix.

    Reference: 
    [1] Chao Chen and Per-Gunnar Martinsson. 2022. Solving linear systems on a GPU 
    with hierarchically off-diagonal low-rank approximations. 
    In Proceedings of the International Conference on High Performance Computing, 
    Networking, Storage and Analysis (SC '22). IEEE Press, Article 84, 1–15.

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format.
    
    b - double
        Array.


    Returns
    --------------------
    x - double
        The solve. 
%}
        
    if isempty(H.D)
        z_alpha = pm_solve(H.A11, b(1: size(H.U1, 1), :));
        z_beta  = pm_solve(H.A22, b(size(H.U1, 1)+1: end, :));

        Y_alpha = pm_solve(H.A11, H.U1);
        Y_beta  = pm_solve(H.A22, H.U2);
        
        z = [z_alpha; z_beta];
        K = blkdiag(H.V1, H.V2) * blkdiag(Y_alpha, Y_beta);

        km1 = size(H.V1,1); km2 = size(Y_alpha, 2);
        K(1:km1, km2+1:end) = eye(km1);
        K(km1+1:end, 1:km2) = eye(km2);

        b = blkdiag(H.V1, H.V2) * z;

        [L, U] = lu(K);
        w = lu_solve(L, U, b);
        x = z - blkdiag(Y_alpha, Y_beta) * w;
    else
        [L, U] = lu(H.D);
        x = lu_solve(L, U, b);
    end
end


function x = lu_solve(L, U, b)
    x = U \ (L \ b);
end