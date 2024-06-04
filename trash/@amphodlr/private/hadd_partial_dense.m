function C = hadd_partial_dense(A, B, operator)
% A must be hodlr format, while B must be dense format
    [mb, nb] = size(B);
    C = zeros(mb, nb);
        
    if isempty(A.D) 
        [m, n, m1, m2, n1, n2] = hsize(A);
        
        if m ~= mb | n ~= nb
            error('Please enter the inputs with consistent dimensions.');
        end  

        C(1:m1, 1:n1) = hadd_partial_dense(A.A11, B(1:m1, 1:n1), operator);
        C(m1+1:end, n1+1:end) = hadd_partial_dense(A.A22, B(m1+1:end, n1+1:end), operator);
        
        if operator == '+'
            C(m1+1:end, 1:n1) = A.U2 * A.V1 + B(m1+1:end, 1:n1);
            C(1:m1, n1+1:end) = A.U1 * A.V2 + B(1:m1, n1+1:end);
        else
            C(m1+1:end, 1:n1) = A.U2 * A.V1 - B(m1+1:end, 1:n1);
            C(1:m1, n1+1:end) = A.U1 * A.V2 - B(1:m1, n1+1:end);
        end

    else
        if operator == '+'
            C = A.D + B;
        else
            C = A.D - B;
        end
    end 
end