function C = hdot_double(A, B)
    if isa(A, 'hodlr') & isa(B, 'hodlr')
        [mA, nA] = hsize(A);
        [mB, nB] = hsize(B);
        if nA ~= mB
            error('Please enter the inputs with correct dimensions.');
        end

        if isempty(A.D) & isempty(B.D)
            A12 = A.U1 * A.V2;
            A21 = A.U2 * A.V1;
            
            B12 = B.U1 * B.V2;
            B21 = B.U2 * B.V1;
                
            C11 = hdot(A.A11, B.A11, 'double') + A12 * B21;
            C12 = hdot(A.A11, B12, 'double') + hdot(A12, B.A22, 'double');
            C21 = hdot(A21, B.A11, 'double') + hdot(A.A22, B21, 'double');
            C22 = A21 * B12 + hdot(A.A22, B.A22, 'double');

            C = [C11, C12; C21, C22];
        
        elseif ~isempty(A.D)
            C = hdot(A.D, B, 'double');
        else
            C = hdot(A, B.D, 'double');
        end
            
    elseif isa(A, 'hodlr')
        if isempty(A.D)
            [mA, nA, su1, su2, sv1, sv2] = hsize(A);
            [mB, nB] = size(B);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.');
            end
        
            C = zeros(mA, nB);
            y1 = hdot(A.A11, B(1:sv1, :), 'double');
            y2 = A.U1 * A.V2 * B(sv1+1:end, :);
            y3 = A.U2 * A.V1 * B(1:sv1, :);
            y4 = hdot(A.A22, B(sv1+1:end, :), 'double');
            
            C(1:su1, :) = y1 + y2;
            C(su1+1:end, :) = y3 + y4;
        else
            C = A.D * B;
        end

    elseif isa(B, 'hodlr')
        if isempty(B.D)
            [mB, nB, su1, su2, sv1, sv2] = hsize(B);
            [mA, nA] = size(A);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.')
            end
        
            C = zeros(mA, nB);
            y1 = hdot(A(:, 1:su1), B.A11, 'double');
            y2 = A(:, su1+1:end) * B.U2 * B.V1;
            y3 = A(:, 1:su1) * B.U1 * B.V2;
            y4 = hdot(A(:, su1+1:end), B.A22, 'double');
            
            C(:, 1:sv1) = y1 + y2;
            C(:, sv1+1:end) = y3 + y4;
        else
            C = A * B.D;
        end
        
    else
        C = A * B;
    end
end
