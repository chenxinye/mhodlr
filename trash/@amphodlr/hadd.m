function C = hadd(varargin)
%{
    The function is calculation of add or minus for HODLR matrix A and B.

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
    
    if (isa(varargin{1}, 'amphodlr') | isa(varargin{1}, 'hodlr') ...
            ) & (isa(varargin{2}, 'amphodlr') | isa(varargin{2}, 'hodlr'))
        input_number = 0;
    elseif isa(varargin{1}, 'hodlr') 
        input_number = 1;
    elseif isa(varargin{2}, 'hodlr') 
        input_number = 2;
    else
        input_number = 3;
    end

    if strcmp(oformat, 'hodlr')
        switch input_number
            case 0
                C = hadd_full_hodlr(varargin{1}, varargin{2}, operator);
            case 1
                C = hadd_partial_hodlr(varargin{1}, varargin{2}, operator);
            case 2
                C = hadd_partial_hodlr(varargin{2}, varargin{1}, operator);
            case 3
                C = hodlr(varargin{1} + varargin{2});
        end
    else
        switch input_number
            case 0
                C = hadd_full_dense(varargin{1}, varargin{2}, operator);
            case 1
                C = hadd_partial_dense(varargin{1}, varargin{2}, operator);
            case 2
                C = hadd_partial_dense(varargin{2}, varargin{1}, operator);
            case 3
                C = varargin{1} + varargin{2};
        end
    end
end