function C = mhdot_dense(A, B)
    if ismember(class(A), {'hodlr', 'amphodlr', 'mphodlr'}) & ismember(class(B), {'hodlr', 'amphodlr', 'mphodlr'})
        [mA, nA] = hsize(A);
        [mB, nB] = hsize(B);
        if nA ~= mB
            error('Please enter the inputs with correct dimensions.');
        end

        if isempty(A.D) & isempty(B.D)
            A12 = mchop(mchop(A.U1) * mchop(A.V2));
            A21 = mchop(mchop(A.U2) * mchop(A.V1));
            
            B12 = mchop(mchop(B.U1) * mchop(B.V2));
            B21 = mchop(mchop(B.U2) * mchop(B.V1));
                
            C11 = mhdot_dense(A.A11, B.A11) + mchop(mchop(A12) * mchop(B21));
            C12 = mhdot_dense(A.A11, B12) + mhdot_dense(A12, B.A22);
            C21 = mhdot_dense(A21, B.A11) + mhdot_dense(A.A22, B21);
            C22 = mchop(mchop(A21) * mchop(B12)) + mhdot_dense(A.A22, B.A22);

            C = [C11, C12; C21, C22];
            
        elseif ~isempty(A.D)
            C = mhdot_dense(A.D, B);
        else
            C = mhdot_dense(A, B.D);
        end
        
            
    elseif ismember(class(A), {'hodlr', 'amphodlr', 'mphodlr'})
        if isempty(A.D)
            [mA, nA, su1, su2, sv1, sv2] = hsize(A);
            [mB, nB] = size(B);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.');
            end
        
            C = zeros(mA, nB);
            y1 = mhdot_dense(A.A11, B(1:sv1, :));
            y2 = mchop(mchop(A.U1) * (mchop(A.V2) * mchop(B(sv1+1:end, :))));
            y3 = mchop(mchop(A.U2) * (mchop(A.V1) * mchop(B(1:sv1, :))));
            y4 = mhdot_dense(A.A22, B(sv1+1:end, :));
            
            C(1:su1, :) = mchop(mchop(y1) + mchop(y2));
            C(su1+1:end, :) = mchop(mchop(y3) + mchop(y4));
        else
            C = mchop(A.D) * mchop(B);
        end

    elseif ismember(class(B), {'hodlr', 'amphodlr', 'mphodlr'})
        if isempty(B.D)
            [mB, nB, su1, su2, sv1, sv2] = hsize(B);
            [mA, nA] = size(A);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.')
            end
        
            C = zeros(mA, nB);
            y1 = mhdot_dense(A(:, 1:su1), B.A11);
            y2 = mchop(mchop(mchop(A(:, su1+1:end)) * mchop(B.U2)) * mchop(B.V1));
            y3 = mchop(mchop(mchop(A(:, 1:su1)) * mchop(B.U1)) * mchop(B.V2));
            y4 = mhdot_dense(A(:, su1+1:end), B.A22);
            
            C(:, 1:sv1) = mchop(y1) + mchop(y2);
            C(:, sv1+1:end) = mchop(y3) + mchop(y4);
        else
            C = mchop(mchop(A) * mchop(B.D));
        end
        
    else
        C = mchop(A) * mchop(B);
    end


    C = mchop(C);
end
