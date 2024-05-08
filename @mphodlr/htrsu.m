function X = htrsu(B, U)
    % Solve XU = B where U is upper triangular matrix
    
    if ~issquare(U)
        error('Please ensure the second input is sqaure matrix.');
    end
    
    if ~(isa(U, 'hodlr') | isa(U, 'mphodlr'))
        error('Please ensure the second input is of a HODLR matrix.');
    end

    if isa(B, 'hodlr') | isa(B, 'mphodlr') 
        % B is of hodlr format
        X = B;
        if isempty(U.D)
            X.A11 = htrsu(B.A11, U.A11);
            X.V1 = htrsu(B.V1, U.A11);

            U12 = U.U1*U.V2;
            [X.U1, X.V2] = compress_m(htrsu(B.U1*B.V2 - hdot(X.A11, U12, 'dense'), U.A22), U.method, U.threshold);
            X.A22 = htrsu(hadd(B.A22, X.U2*X.V1*U12, '-'), U.A22);
        else
            X.D = mrdivide(B.D, U.D);
        end
        
    else
        % B is dense format 
        if isempty(U.D)
            
            m1 = hsize(U, 2);
            X11 = htrsu(B(1:m1, 1:m1), U.A11);
            X21 = htrsu(B(m1+1:end, 1:m1), U.A11);
            
            U12 = U.U1*U.V2;
            X12 = htrsu(B(1:m1, m1+1:end) - X11 * U12, U.A22);
            X22 = htrsu(B(m1+1:end, m1+1:end) - X21*U12, U.A22);
            X = [X11, X12; X21, X22];
        else
            X = mrdivide(B, U.D);
        end
    end


end