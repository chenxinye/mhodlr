function X = mhtrsu(B, U, varargin)
    %{
        Parameters
        --------------------
        B - hodlr | double | single
            Matrix in HODLR format - hodlr class.
    
        U - hodlr | double | single
            Matrix in HODLR format - hodlr class.
    
        itype- int
            `1`: Solve XU = B, U is upper triangular matrix
            
            `2` Solve BX = U, B is upper triangular matrix.
    
    
        Returns
        --------------------
        X - double 
            The solve.
    
    %}
    
    
        if nargin == 3
            itype = varargin{1}; % if itype is 1, then the method is used for LU factorization
        else 
            itype = 1;
        end
    
        if itype ~= 2 % Solve XU = B where U is upper triangular matrix
            if ~issquare(U)
                error('Please ensure the second input is sqaure matrix.');
            end
            
            if ~(isa(U, 'hodlr') | isa(U, 'mphodlr') | isa(U, 'amphodlr'))
                error('Please ensure the second input is of a HODLR matrix.');
            end
            
            U = hmchop(U);

            if isa(B, 'hodlr') | isa(B, 'mphodlr') | isa(B, 'amphodlr') 
                % B is of hodlr format
                B = hmchop(B);
                X = B;
                if isempty(U.D)
                    X.A11 = mhtrsu(B.A11, U.A11);
                    X.V1 = mhtrsu(B.V1, U.A11, 1);
    
                    U12 = mchop(U.U1*U.V2);
                    [X.U1, X.V2, ~] = compress_m(mhtrsu(B.U1*B.V2 - hdot(X.A11, U12, 'dense'), U.A22, 1), U.method, U.vareps, U.max_rnk, U.trun_norm_tp, U.issparse);
                    
                    X.U1 = mchop(X.U1);
                    X.V2 = mchop(X.V2);
                    X.A22 = mhtrsu(mchop(hadd(B.A22, mchop(X.U2*mchop(X.V1*U12)), '-')), U.A22);
                else
                    X.D = mchop(mrdivide(B.D, U.D));
                end
                
            else
                % B is dense format 
                B = mchop(B);
                if isempty(U.D) % XU = B
                    [~, ~, n1] = size_t(U, 2);
                    X1 = mhtrsu(B(:, 1:n1), U.A11);
                    U12 = mchop(U.U1*U.V2);
                    X2 = mhtrsu(B(:, n1+1:end) - X1 * U12, U.A22);
                    X = mchop([X1, X2]);
                else
                    X = mchop(mrdivide(B, U.D));
                end
            end
        else
            % BX = U
            if isempty(B.D)
                m1 = size_t(B, 2);
                X2 = mhtrsu(B.A22, U(m1+1:end, :), 2);
                B12 = mchop(B.U1*B.V2);
                X1 = mhtrsu(B.A11, U(1:m1, :) - B12*X2, 2);
                X = mchop([X1; X2]);
            else
                X = mchop(mldivide(B.D, U));
            end
        end
    end