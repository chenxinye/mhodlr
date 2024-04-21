function C = hdot(A, B) % one if left A, the other is right A
    if isa(A, 'hodlr') & isa(B, 'hodlr')
        error('Please enter correct types of inputs.');
    end
    
    if isa(A, 'hodlr')
        if isempty(A.D)
            su1 = size(A.U1, 1);
            su2 = size(A.U2, 1);
            sv1 = size(A.V1, 2);
            sv2 = size(A.V2, 2);
        
            mHRows = su1 + su2;
            nHCols = sv1 + sv2;
            
            [mRows, nCols] = size(B);
            
            if nHCols ~= mRows
                error('Please enter the inputs with correct dimensions.');
            end
        
            C = zeros(mHRows, nCols);
            y1 = hdot(A.A11, B(1:su1, :));
            y2 = A.U1 * A.V2 * B(su1+1:end, :);
            y3 = A.U2 * A.V1 * B(1:su1, :);
            y4 = hdot(A.A22, B(su1+1:end, :));
            
        
            C(1:su1, :) = y1 + y2;
            C(su1+1:end, :) = y3 + y4;
        else
            C = A.D * B;
        end

    elseif isa(B, 'hodlr')
        if isempty(B.D)
            su1 = size(B.U1, 1);
            su2 = size(B.U2, 1);
            sv1 = size(B.V1, 2);
            sv2 = size(B.V2, 2);
        
            mHRows = su1 + su2;
            nHCols = sv1 + sv2;
            
            [mRows, nCols] = size(A);
            
            if nCols ~= mHRows
                error('Please enter the inputs with correct dimensions.')
            end
        
            C = zeros(mRows, nHCols);
            y1 = hdot(A(:, 1:su1), B.A11);
            y2 = A(:, su1+1:end) * B.U2 * B.V1;
            y3 = A(:, 1:su1) * B.U1 * B.V2;
            y4 = hdot(A(:, su1+1:end), B.A22);
            
            C(:, 1:su1) = y1 + y2;
            C(:, su1+1:end) = y3 + y4;
        else
            C = A * B.D;
        end
        
    else
        C = A * B;
    end
end