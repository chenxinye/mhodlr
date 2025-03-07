function C = hadd_full_dense(A, B, operator)
% Output dense format when inputs A, B are HODLR matrix.
% 
    if isempty(A.D) & isempty(B.D) 
        [m, n] = size_t(A);
        m1 = size(A.U1, 1);
        n1 = size(A.V1, 2);

        [mb, nb] = size_t(B);
        C = zeros(mb, nb);

        if m ~= mb | n ~= nb
            error('Please enter the inputs with consistent dimensions.');
        end  

        C(1:m1, 1:n1) = hadd_full_dense(A.A11, B.A11, operator);
        C(m1+1:end, n1+1:end) = hadd_full_dense(A.A22, B.A22, operator);
        
        if operator == '+'
            C(m1+1:end, 1:n1) = A.U2 * A.V1 + B.U2 * B.V1;
            C(1:m1, n1+1:end) = A.U1 * A.V2 + B.U1 * B.V2;
        else
            C(m1+1:end, 1:n1) = A.U2 * A.V1 - B.U2 * B.V1;
            C(1:m1, n1+1:end) = A.U1 * A.V2 - B.U1 * B.V2;
        end

    else
        if operator == '+'
            C = A.D + B.D;
        else
            C = A.D - B.D;
        end
    end 
end