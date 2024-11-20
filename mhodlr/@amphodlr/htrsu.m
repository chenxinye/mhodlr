function X = htrsu(B, U, varargin)
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
    
        if isa(B, 'hodlr') | isa(B, 'mphodlr')  | isa(B, 'amphodlr') 
            % B is of hodlr format
            X = B;
            if isempty(U.D)
                X.A11 = htrsu(B.A11, U.A11);
                X.V1 = htrsu(B.V1, U.A11, 1);

                U12 = U.U1*U.V2;
                [X.U1, X.V2] = compress_m(htrsu(B.U1*B.V2 - hdot(X.A11, U12, 'dense'), U.A22, 1), U.method, U.vareps, obj.issparse);
            
                X.A22 = htrsu(hadd(B.A22, X.U2*X.V1*U12, '-'), U.A22);
            else
                X.D = mrdivide(B.D, U.D);
            end
            
        else
            % B is dense format 
            if isempty(U.D) % XU = B
                [~, ~, n1] = hsize(U, 2);
                X1 = htrsu(B(:, 1:n1), U.A11);
                U12 = U.U1*U.V2;
                X2 = htrsu(B(:, n1+1:end) - X1 * U12, U.A22);
                X = [X1, X2];

                % [m, n] = size(B);
                % if itype == 1
                %     X11 = htrsu(B(1:m1, 1:n1), U.A11, 1);
                %     X21 = htrsu(B(m1+1:end, 1:m1), U.A11, 1);
                %     
                %     U12 = U.U1*U.V2;
                %     X12 = htrsu(B(1:m1, m1+1:end) - X11 * U12, U.A22, 1);
                %     X22 = htrsu(B(m1+1:end, m1+1:end) - X21 * U12, U.A22, 1);
                %     X = [X11, X12; X21, X22];
                % else
                %     X1 = htrsu(B(:, 1:m1), U.A11);
                %     U12 = U.U1*U.V2;
                %     X2 = htrsu(B(:, m1+1:end) - X1 * U12, U.A22);
                %     X = [X1, X2];
                % end
            else
                X = mrdivide(B, U.D);
            end
        end
    else
        % BX = U
        if isempty(B.D)
            m1 = hsize(B, 2);
            X2 = htrsu(B.A22, U(m1+1:end, :), 2);
            B12 = B.U1*B.V2;
            X1 = htrsu(B.A11, U(1:m1, :) - B12*X2, 2);
            X = [X1; X2];
        else
            X = mldivide(B.D, U);
        end
    end
end
