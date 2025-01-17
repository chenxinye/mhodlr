function C = hdot_hodlr(A, B)
    if A.bottom_level ~= B.bottom_level
        error("Please ensure the two input matrices are of the same cluster structure.")
    else
        vareps = max(A.vareps, B.vareps);
    end

    if ~isempty(A.D)
        D = A.D * B.D;
        C = hodlr(D, 0, A.min_block_size);
        C.vareps = vareps;
        [m, n] = size(D);
        C.shape = [m, n];
    else
        C = A;
        C.vareps = vareps;
        C.A11 = hadd_partial_hodlr(hdot_hodlr(A.A11, B.A11), A.U1*(A.V2*B.U2)*B.V1, '+', true);
        C.A22 = hadd_partial_hodlr(hdot_hodlr(A.A22, B.A22), A.U2*(A.V1*B.U1)*B.V2, '+', true);

        
        A12 = hdot_dense(A.A11, B.U1*B.V2) + hdot_dense(A.U1*A.V2, B.A22);
        A21 = hdot_dense(A.U2*A.V1, B.A11) + hdot_dense(A.A22, B.U2*B.V1);
        [C.U1, C.V2] = compress_m(A12, 'svd', vareps);
        [C.U2, C.V1] = compress_m(A21, 'svd', vareps);
    end
end