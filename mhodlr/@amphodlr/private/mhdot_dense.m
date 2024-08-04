function C = mhdot_dense(A, B)
    if ismember(class(A), {'hodlr', 'amphodlr', 'mphodlr'}) & ismember(class(B), {'hodlr', 'amphodlr', 'mphodlr'})
        [mA, nA] = hsize(A);
        [mB, nB] = hsize(B);

        A = hmchop(A);
        B = hmchop(B);

        if nA ~= mB
            error('Please enter the inputs with correct dimensions.');
        end

        if isempty(A.D) & isempty(B.D)
            A12 = mchop(A.U1 * A.V2);
            A21 = mchop(A.U2 * A.V1);
            
            B12 = mchop(B.U1 * B.V2);
            B21 = mchop(B.U2 * B.V1);
                
            C11 = mhdot_dense(A.A11, B.A11) + mchop(A12 * B21);
            C12 = mhdot_dense(A.A11, B12) + mhdot_dense(A12, B.A22);
            C21 = mhdot_dense(A21, B.A11) + mhdot_dense(A.A22, B21);
            C22 = mchop(A21 * B12) + mhdot_dense(A.A22, B.A22);

            C = [C11, C12; C21, C22];
            
        elseif ~isempty(A.D)
            C = mhdot_dense(A.D, B);
        else
            C = mhdot_dense(A, B.D);
        end
        
            
    elseif ismember(class(A), {'hodlr', 'amphodlr', 'mphodlr'})
        A = hmchop(A);
        B = mchop(B);

        if isempty(A.D)
            [mA, nA, su1, ~, sv1, ~] = hsize(A);
            [mB, nB] = size(B);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.');
            end
        
            C = zeros(mA, nB);
            y1 = mhdot_dense(A.A11, B(1:sv1, :));
            y2 = mchop(A.U1 * mchop(A.V2 * B(sv1+1:end, :)));
            y3 = mchop(A.U2 * mchop(A.V1 * B(1:sv1, :)));
            y4 = mhdot_dense(A.A22, B(sv1+1:end, :));
            
            C(1:su1, :) = mchop(y1 + y2);
            C(su1+1:end, :) = mchop(y3 + y4);


        else
            C = mchop(mchop(A.D) * mchop(B));
        end

    elseif ismember(class(B), {'hodlr', 'amphodlr', 'mphodlr'})
        A = mchop(A);
        B = hmchop(B);

        if isempty(B.D)
            [mB, nB, su1, ~, sv1, ~] = hsize(B);
            [mA, nA] = size(A);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.')
            end
        
            C = zeros(mA, nB);
            y1 = mhdot_dense(A(:, 1:su1), B.A11);
            y2 = mchop(mchop(A(:, su1+1:end) * B.U2) * B.V1);
            y3 = mchop(mchop(A(:, 1:su1) * B.U1) * B.V2);
            y4 = mhdot_dense(A(:, su1+1:end), B.A22);
            
            C(:, 1:sv1) = mchop(y1 + y2);
            C(:, sv1+1:end) = mchop(y3 + y4);
        else
            C = mchop(A * B.D);
        end
        
    else
        C = mchop(mchop(A) * mchop(B));
    end

end
