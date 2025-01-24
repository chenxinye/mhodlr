function C = mhadd_full_dense(A, B, operator)
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

        C(1:m1, 1:n1) = mchop(mhadd_full_dense(A.A11, B.A11, operator));
        C(m1+1:end, n1+1:end) = mchop(mhadd_full_dense(A.A22, B.A22, operator));
        
        if operator == '+'
            C(m1+1:end, 1:n1) = mchop(mchop(A.U2 * A.V1) + mchop(B.U2 * B.V1));
            C(1:m1, n1+1:end) = mchop(mchop(A.U1 * A.V2) + mchop(B.U1 * B.V2));
        else
            C(m1+1:end, 1:n1) = mchop(mchop(A.U2 * A.V1) - mchop(B.U2 * B.V1));
            C(1:m1, n1+1:end) = mchop(mchop(A.U1 * A.V2) - mchop(B.U1 * B.V2));
        end

    else
        if operator == '+'
            C = mchop(A.dense + B.dense);
        else
            C = mchop(A.dense - B.dense);
        end
    end 
end