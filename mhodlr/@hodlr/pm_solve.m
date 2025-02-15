function x = pm_solve(H, b)

    if isempty(H.D)
        z_alpha = pm_solve(H.A11, b(1: size(H.U1, 1)));
        z_beta  = pm_solve(H.A22, b(size(H.U1, 1)+1: end));

        Y_alpha = pm_solve(H.A11, H.U1);
        Y_beta  = pm_solve(H.A22, H.U2);
        
        z = [z_alpha; z_beta];
        K = blkdiag(H.V1, H.V2) * blkdiag(Y_alpha, Y_beta);

        K(1:size(H.V1,1), size(Y_alpha,2)+1:end) = eye(size(H.V1, 1));
        K(size(H.V1, 1)+1:end, 1:size(Y_alpha,2)) = eye(size(Y_beta, 2));

        [L, U] = lu(K);
        b = blkdiag(H.V1, H.V2) * z;
        w = lu_solve(L, U, b);
    
        x = z - blkdiag(Y_alpha, Y_beta) * w;
        
    else
        x = H.D\b;
    end

end


function x = lu_solve(L, U, b)

    x = U \ L  \ b ;
end