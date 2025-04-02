function X = mhtrsu(B, U, varargin)
%{
    Parameters
    --------------------
    B - hodlr | double | single
        Matrix in HODLR format - hodlr class.

    U - hodlr | double | single
        Matrix in HODLR format - hodlr class.

    itype- int
        `1`: Solve XU = B, U is upper triangular matrix, 
              e.g., QR = A
        
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
    % Compute class flags once
    isU_hodlr = is_hodlr_class(U);
    isB_hodlr = is_hodlr_class(B);

    if itype ~= 2 % Solve XU = B where U is upper triangular matrix
        if ~issquare(U)
            error('Please ensure the second input is sqaure matrix.');
        end
        
        if ~isU_hodlr
            error('Please ensure the second input is of a HODLR matrix.');
        end
    
        U = hmchop(U);

        if isB_hodlr
            B = hmchop(B);
            % B is of hodlr format
            X = B;
            if isempty(U.D)
                X.A11 = mhtrsu(B.A11, U.A11);
                X.V1 = mhtrsu(B.V1, U.A11, 1);

                [X.U1, X.V2, ~] = mp_compress_m(mhtrsu(B.U1*B.V2 - mhdot(X.A11, U.U1, 'dense') * U.V2, U.A22, 1), U.method, U.vareps, U.max_rnk, U.trun_norm_tp, U.issparse);
            
                X.A22 = mhtrsu(mhadd(B.A22, mchop(mchop(X.U2*mchop(X.V1*U.U1))*U.V2), '-'), U.A22);
            else
                X.D = mchop(mrdivide(B.dense, U.D));
            end
            
        else
            % B is dense format 
            B = mchop(B);
            if isempty(U.D) % XU = B
                [~, ~, n1] = size_t(U, 2);
                X1 = mhtrsu(B(:, 1:n1), U.A11);
                X2 = mhtrsu(B(:, n1+1:end) - mchop(mchop(X1 * U.U1) * U.V2), U.A22);
                X = [X1, X2];

                % [m, n] = size(B);
                % if itype == 1
                %     X11 = mhtrsu(B(1:m1, 1:n1), U.A11, 1);
                %     X21 = mhtrsu(B(m1+1:end, 1:m1), U.A11, 1);
                %     
                %     U12 = U.U1*U.V2;
                %     X12 = mhtrsu(B(1:m1, m1+1:end) - X11 * U12, U.A22, 1);
                %     X22 = mhtrsu(B(m1+1:end, m1+1:end) - X21 * U12, U.A22, 1);
                %     X = [X11, X12; X21, X22];
                % else
                %     X1 = mhtrsu(B(:, 1:m1), U.A11);
                %     U12 = U.U1*U.V2;
                %     X2 = mhtrsu(B(:, m1+1:end) - X1 * U12, U.A22);
                %     X = [X1, X2];
                % end
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
            X1 = mhtrsu(B.A11, mchop(U(1:m1, :) - mchop(B12*X2)), 2);
            X = [X1; X2];
        else
            X = mchop(mldivide(B.D, U));
        end
    end
end


function is_hodlr = is_hodlr_class(obj)
    % Helper function to check if obj is a HODLR-like class
    switch class(obj)
        case {'hodlr', 'amphodlr', 'mphodlr'}
            is_hodlr = true;
        otherwise
            is_hodlr = false;
    end
end
