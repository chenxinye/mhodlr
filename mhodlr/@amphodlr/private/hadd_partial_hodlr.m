function C = hadd_partial_hodlr(A, B, operator, B_is_dense)
    % A must be hodlr format, while B must be dense format

    if B_is_dense
        [mB, nB] = size(B);
        C = A;

        if isempty(A.D) 
            [m, n, m1, ~, n1, ~] = size_t(A);
            
            if mB == 0 | nB == 0 | m == 0 | n == 0
                C = hodlr;
            else
                if m ~= mB | n ~= nB
                    error('Please enter the inputs with consistent dimensions.');
                end  
                
                C.A11 = hadd_partial_hodlr(A.A11, B(1:m1, 1:n1), operator, B_is_dense);
                C.A22 = hadd_partial_hodlr(A.A22, B(m1+1:end, n1+1:end), operator, B_is_dense);

                if operator == '+'
                    [C.U2, C.V1, ~] = A.compress(A.U2 * A.V1 + B(m1+1:end, 1:n1));
                    [C.U1, C.V2, ~] = A.compress(A.U1 * A.V2 + B(1:m1, n1+1:end));
                else
                    [C.U2, C.V1, ~] = A.compress(A.U2 * A.V1 - B(m1+1:end, 1:n1));
                    [C.U1, C.V2, ~] = A.compress(A.U1 * A.V2 - B(1:m1, n1+1:end));
                end
            end
        else
            if operator == '+'
                C.D = A.D + B;
            else
                C.D = A.D - B;
            end
        end
    else
        C = B;
        if isempty(B.D) 
            [m, n] = size(A);
            [mB, nB, m1, ~, n1, ~] = size_t(B);

            if mB == 0 | nB == 0 | m == 0 | n == 0
                C = hodlr;
            else
                if m ~= mB | n ~= nB
                    error('Please enter the inputs with consistent dimensions.');
                end  
                
                C.A11 = hadd_partial_hodlr(A(1:m1, 1:n1), B.A11, operator, B_is_dense);
                C.A22 = hadd_partial_hodlr(A(m1+1:end, n1+1:end), B.A22, operator, B_is_dense);

                if operator == '+'
                    [C.U2, C.V1, ~] = B.compress(A(m1+1:end, 1:n1) + B.U2 * B.V1);
                    [C.U1, C.V2, ~] = B.compress(A(1:m1, n1+1:end) + B.U1 * B.V2);
                else
                    [C.U2, C.V1, ~] = B.compress(A(m1+1:end, 1:n1) - B.U2 * B.V1);
                    [C.U1, C.V2, ~] = B.compress(A(1:m1, n1+1:end) - B.U1 * B.V2);
                end

            end
        else
            if operator == '+'
                C.D = A + B.dense;
            else
                C.D = A- B.dense;
            end
        end
    end
end