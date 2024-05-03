function C = mphdot_double(A, B)
    if isa(A, 'hodlr') & isa(B, 'hodlr')
        [mA, nA] = hsize(A);
        [mB, nB] = hsize(B);
        if nA ~= mB
            error('Please enter the inputs with correct dimensions.');
        end

        set_prec(A.prec_settings{A.level});
        
        if isempty(A.D) & isempty(B.D)
            A12 = mchop(mchop(A.U1) * mchop(A.V2));
            A21 = mchop(mchop(A.U2) * mchop(A.V1));
            
            B12 = mchop(mchop(B.U1) * mchop(B.V2));
            B21 = mchop(mchop(B.U2) * mchop(B.V1));
                
            C11 = hdot_double(A.A11, B.A11, 'double') + mchop(A12 * B21);
            C12 = hdot_double(A.A11, B12, 'double') + hdot_double(A12, B.A22, 'double');
            C21 = hdot_double(A21, B.A11, 'double') + hdot_double(A.A22, B21, 'double');
            C22 = mchop(A21 * B12) + hdot_double(A.A22, B.A22, 'double');

            C = [C11, C12; C21, C22];
        
        elseif ~isempty(A.D)
            C = hdot_double(A.D, B, 'double');
        else
            C = hdot_double(A, B.D, 'double');
        end
            
    elseif isa(A, 'hodlr')
        set_prec(A.prec_settings{A.level});

        if isempty(A.D)
            [mA, nA, su1, su2, sv1, sv2] = hsize(A);
            [mB, nB] = size(B);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.');
            end
        
            C = zeros(mA, nB);
            y1 = hdot_double(A.A11, B(1:sv1, :), 'double');
            y2 = mchop(mchop(mchop(A.U1) * mchop(A.V2)) * mchop(B(sv1+1:end, :)));
            y3 = mchop(mchop(mchop(A.U2) * mchop(A.V1)) * mchop(B(1:sv1, :)));
            y4 = hdot_double(A.A22, B(sv1+1:end, :), 'double');
            
            C(1:su1, :) = mchop(y1 + y2);
            C(su1+1:end, :) = mchop(y3 + y4);
        else
            C = mchop(mchop(A.D) * mchop(B));
        end

    elseif isa(B, 'hodlr')
        set_prec(B.prec_settings{B.level});

        if isempty(B.D)
            [mB, nB, su1, su2, sv1, sv2] = hsize(B);
            [mA, nA] = size(A);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.')
            end
        
            C = zeros(mA, nB);
            y1 = hdot_double(A(:, 1:su1), B.A11, 'double');
            y2 = mchop(mchop(mchop(A(:, su1+1:end)) * mchop(B.U2)) * mchop(B.V1));
            y3 = mchop(mchop(mchop(A(:, 1:su1)) * B.U1) * B.V2);
            y4 = hdot_double(A(:, su1+1:end), B.A22, 'double');
            
            C(:, 1:sv1) = mchop(y1 + y2);
            C(:, sv1+1:end) = mchop(y3 + y4);
        else
            C = mchop(mchop(A) * mchop(B.D));
        end
        
    else
        C = mchop(mchop(A) * mchop(B));
    end
end
