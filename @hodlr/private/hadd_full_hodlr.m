function C = hadd_full_hodlr(A, B, operator)
    [m, n] = hsize(A);
    [mb, nb] = hsize(B);
    C = A;

    if isempty(A.D) & isempty(B.D) 
        if m ~= mb | n ~= nb
            error('Please enter the inputs with consistent dimensions.');
        end  

        C.A11 = hadd_full_hodlr(A.A11, B.A11, operator);
        C.A22 = hadd_full_hodlr(A.A22, B.A22, operator);
        
        if operator == '+'
            [C.U2, C.V1] = A.compress(A.U2 * A.V1 + B.U2 * B.V1);
            [C.U1, C.V2] = A.compress(A.U1 * A.V2 + B.U1 * B.V2);
        else
            [C.U2, C.V1] = A.compress(A.U2 * A.V1 - B.U2 * B.V1);
            [C.U1, C.V2] = A.compress(A.U1 * A.V2 - B.U1 * B.V2);
        end

    else
        if operator == '+'
            C.D = A.D + B.D;
        else
            C.D = A.D - B.D;
        end
    end 
end