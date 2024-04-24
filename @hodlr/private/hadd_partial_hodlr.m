function C = hadd_partial_hodlr(A, B, operator)
% A must be hodlr format, while B must be double format
    [mb, nb] = size(B);
    C = A;

    if isempty(A.D) 
        [m, n, m1, m2, n1, n2] = hsize(A);
        
        if m ~= mb | n ~= nb
            error('Please enter the inputs with consistent dimensions.');
        end  
        
        C.A11 = hadd_partial_hodlr(A.A11, B(1:m1, 1:n1), operator);
        C.A22 = hadd_partial_hodlr(A.A22, B(m1+1:end, n1+1:end), operator);

        if operator == '+'
            [C.U2, C.V1] = A.compress(A.U2 * A.V1 + B(m1+1:end, 1:n1));
            [C.U1, C.V2] = A.compress(A.U1 * A.V2 + B(1:m1, n1+1:end));
        else
            [C.U2, C.V1] = A.compress(A.U2 * A.V1 - B(m1+1:end, 1:n1));
            [C.U1, C.V2] = A.compress(A.U1 * A.V2 - B(1:m1, n1+1:end));
        end
    else
        if operator == '+'
            C.D = A.D + B;
        else
            C.D = A.D - B;
        end
    end
end