function x = pm_solve(H, b)

    if ~isempty(H.D)
        midm = ceil(H.shape[0] / 2)
        z_alpha = pm_solve(H.A11, b_alpha(1: midm));
        z_beta  = pm_solve(H.A22, b_beta(midm+1: end));

        Y_alpha = pm_solve(H.A11, H.U1);
        Y_beta  = pm_solve(H.A22, H.U2);
        z = [z_alpha; z_beta];

        K = blkdiag(H.V1 * Y_alpha, H.V2 * Y_beta);
        
        K(1:midm, midn+1:) = eye(midm);
        K(midm+1:, 1:midm) = eye(midm);

        L, U = lu(K);
        b = blkdiag(H.V1, H.V2) * z;
        w = lu_solve(L, U, b);
    
        x = z - blkdiag(Y_alpha, Y_beta) * w;
        
    else
        L, U = lu(H.D);
        x = lu_solve(L, U, b);
    end

end


function x = lu_solve(L, U, b)
    x = b / L / U;
end