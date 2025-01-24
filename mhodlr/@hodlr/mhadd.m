function C = mhadd(varargin)
%{
    The function is used for the operation of summation or subtraction for HODLR matrix A and B.

    Parameters
    --------------------
    A - hodlr | double
        Input matrix - hodlr class / dense tyle.
    
    B - hodlr | double
        Input matrix - hodlr class / dense tyle.
    
    operator - str, default = '+'
        The operator of add ('+') or minus ('-'), string type.

    oformat - str, default = 'hodlr'
        The format of returns.
    
    Returns
    --------------------
    C - hodlr | double
        Return matrix in hodlr class or dense array.
    
%}

    if nargin == 2
        operator = '+';
        oformat = 'hodlr';

    elseif nargin == 3
        operator = varargin{3};
        oformat = 'hodlr';

    elseif nargin == 4 
        operator = varargin{3};
        oformat = varargin{4};

    elseif nargin > 4
        error('Please enter the correct number of inputs.');
    end


    A = varargin{1};
    B = varargin{2};
    
    if ismember(class(A), {'hodlr', 'amphodlr', 'mphodlr'}) & ismember(class(B), {'hodlr', 'amphodlr', 'mphodlr'})
        input_number = 0;
        A = hmchop(A);
        B = hmchop(B);
    elseif ismember(class(A), {'hodlr', 'amphodlr', 'mphodlr'})
        input_number = 1;
        A = hmchop(A);
    elseif ismember(class(B), {'hodlr', 'amphodlr', 'mphodlr'})
        input_number = 2;
        B = hmchop(B);
    else
        input_number = 3;
        A = mchop(A);
        B = mchop(B);
    end

    if strcmp(oformat, 'hodlr')
        switch input_number
            case 0
                C = mhadd_full_hodlr(A, B, operator);
            case 1
                C = mhadd_partial_hodlr(A, B, operator, true);
            case 2
                C = mhadd_partial_hodlr(A, B, operator, false);
            case 3
                if strcmp(operator, '-')
                    C = hodlr(A + B);
                else
                    C = hodlr(A - B);
                end

                C = hmchop(C);
        end

    else
        switch input_number
            case 0
                C = mhadd_full_dense(A, B, operator);
            case 1
                C = mhadd_partial_dense(A, B, operator);
            case 2
                C = mhadd_partial_dense(B, A, operator);
                
                if strcmp(operator, '-')
                    C = -C;
                end
            case 3
                if strcmp(operator, '-')
                    C = A + B;
                else
                    C = A - B;
                end
        end
    end
end