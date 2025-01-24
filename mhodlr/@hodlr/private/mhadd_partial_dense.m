function C = mhadd_partial_dense(A, B, operator)
% A must be hodlr format, while B must be dense format
    [mb, nb] = size(B);
    C = zeros(mb, nb);
        
    if isempty(A.D) 
        [m, n, m1, ~, n1, ~] = size_t(A);
        
        if m ~= mb | n ~= nb
            error('Please enter the inputs with consistent dimensions.');
        end  

        C(1:m1, 1:n1) = mchop(mhadd_partial_dense(A.A11, B(1:m1, 1:n1), operator));
        C(m1+1:end, n1+1:end) = mchop(mhadd_partial_dense(A.A22, B(m1+1:end, n1+1:end), operator));
        
        if operator == '+'
            C(m1+1:end, 1:n1) = mchop(mchop(A.U2 * A.V1) + mchop(B(m1+1:end, 1:n1)));
            C(1:m1, n1+1:end) = mchop(mchop(A.U1 * A.V2) + mchop(B(1:m1, n1+1:end)));
        else
            C(m1+1:end, 1:n1) = mchop(mchop(A.U2 * A.V1) - mchop(B(m1+1:end, 1:n1)));
            C(1:m1, n1+1:end) = mchop(mchop(A.U1 * A.V2) - mchop(B(1:m1, n1+1:end)));
        end

    else
        if operator == '+'
            C = mchop(A.dense + B);
        else
            C = mchop(A.dense - B);
        end
    end 
end