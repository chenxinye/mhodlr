function X = htrsl(L, B, varargin)
    % Solve LX = B where L is lower triangular matrix

    if ~issquare(L)
        error('Please ensure the first input is sqaure matrix.');
    end
    
    if ~(isa(L, 'hodlr') | isa(L, 'mphodlr'))
        error('Please ensure the first input is of a HODLR matrix.');
    end

    if nargin > 2
        itype = varargin{1}; % if itype is 1, then the method is used for LU factorization
    else
        itype = 0;
    end

    if isa(B, 'hodlr') | isa(B, 'mphodlr') 
        % B is of hodlr format
        X = B;
        if isempty(L.D)
            X.A11 = htrsl(L.A11, B.A11);
            X.U1 = htrsl(L.A11, B.U1, 1);

            B21 = B.U2*B.V1;
            L21 = L.U2*L.V1;
            [X.U2, X.V1] = compress_m(htrsl(L.A22, B21 - hdot(L21, X.A11, 'dense')), L.method, L.threshold);
            X.A22 = htrsl(L.A22, hadd(B.A22, L21*X.U1*X.V2, '-'));
        else
            X.D = mldivide(L.D, B.D);
        end
    
    else
        % B is dense format 
        if isempty(L.D) % think B is a block matrix
            m1 = hsize(L, 2);
            if itype == 1
                [m, n] = hsize(L, 1);
                [a, b] = size(B);
                
                X11 = htrsl(L.A11, B(1:m1, 1:m1));
                X12 = htrsl(L.A11, B(1:m1, m1+1:end));
                L21 = L.U2*L.V1;

                X21 = htrsl(L.A22, B(m1+1:end, 1:m1) - L21 * X11);
                X22 = htrsl(L.A22, B(m1+1:end, m1+1:end) - L21 * X12);
                X = [X11, X12; X21, X22];
            else
                X1 = htrsl(L.A11, B(1:m1, :));
                L21 = L.U2*L.V1;
                X2 = htrsl(L.A22, B(m1+1:end, :) - L21 * X1);
                X = [X1; X2];
            end
        else
            X = mldivide(L.D, B);
        end
    end


end