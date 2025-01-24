function C = hadd(varargin)
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
    
    if (isa(A, 'hodlr') | isa(A, 'mphodlr') | isa(A, 'amphodlr') ...
            ) & (isa(B, 'hodlr') | isa(B, 'mphodlr') | isa(B, 'amphodlr'))
        input_number = 0;
    elseif isa(A, 'hodlr') | isa(A, 'mphodlr') | isa(A, 'amphodlr') 
        input_number = 1;
    elseif isa(B, 'hodlr') | isa(B, 'mphodlr') | isa(B, 'amphodlr') 
        input_number = 2;
    else
        input_number = 3;
    end

    if strcmp(oformat, 'hodlr')
        switch input_number
            case 0
                C = hadd_full_hodlr(A, B, operator);
            case 1
                C = hadd_partial_hodlr(A, B, operator, true);
            case 2
                C = hadd_partial_hodlr(A, B, operator, false);
            case 3
                if strcmp(operator, '-')
                    C = hodlr(A + B);
                else
                    C = hodlr(A - B);
                end
        end
    else
        switch input_number
            case 0
                C = hadd_full_dense(A, B, operator);
            case 1
                C = hadd_partial_dense(A, B, operator);
            case 2
                C = hadd_partial_dense(B, A, operator);
                
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