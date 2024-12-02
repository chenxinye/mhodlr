function X = mhtrsl(L, B, varargin)
    % Solve LX = B where L is lower triangular matrix

    if ~issquare(L)
        error('Please ensure the first input is sqaure matrix.');
    end
    
    if ~(isa(L, 'hodlr') | isa(L, 'mphodlr') | isa(L, 'amphodlr'))
        error('Please ensure the first input is of a HODLR matrix.');
    end

    L = hmchop(L);

    if nargin > 2
        itype = varargin{1}; % if itype is 1, then the method is used for LU factorization
    else
        itype = 0;
    end

    if isa(B, 'hodlr') | isa(B, 'mphodlr') | isa(B, 'amphodlr')
        B = hmchop(B);
        % B is of hodlr format
        X = B;
        if isempty(L.D)
            X.A11 = mchop(mhtrsl(L.A11, B.A11));
            X.U1 = mchop(mhtrsl(L.A11, B.U1, 1));

            B21 = mchop(B.U2*B.V1);
            L21 = mchop(L.U2*L.V1);
            [X.U2, X.V1] = compress_m(mhtrsl(L.A22, B21 - hdot(L21, X.A11, 'dense')), L.method, L.vareps, L.max_rnk, L.trun_norm_tp, L.issparse);
            X.A22 = mhtrsl(L.A22, hadd(B.A22, L21*X.U1*X.V2, '-'));
            X = hmchop(X);
        else
            X.D = mchop(mldivide(L.D, B.D));
        end
    
    else
        B = mchop(B);
        % B is dense format 
        if isempty(L.D) % think B is a block matrix
            m1 = hsize(L, 2);
            if itype == 1
                [m, n] = hsize(L, 1);
                [a, b] = size(B);
                
                X11 = mchop(mhtrsl(L.A11, B(1:m1, 1:m1)));
                X12 = mchop(mhtrsl(L.A11, B(1:m1, m1+1:end)));
                L21 = mchop(L.U2*L.V1);

                X21 = mchop(mhtrsl(L.A22, mchop(B(m1+1:end, 1:m1) - mchop(L21 * X11))));
                X22 = mchop(mhtrsl(L.A22, mchop(B(m1+1:end, m1+1:end) - mchop(L21 * X12))));
                X = [X11, X12; X21, X22];
            else
                X1 = mchop(mhtrsl(L.A11, B(1:m1, :)));
                L21 = mchop(L.U2*L.V1);
                X2 = mchop(mhtrsl(L.A22, mchop(B(m1+1:end, :) - mchop(L21 * X1))));
                X = [X1; X2];
            end
        else
            X = mchop(mldivide(L.D, B));
        end
    end


end
