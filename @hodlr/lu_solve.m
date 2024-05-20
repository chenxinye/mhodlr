function x = lu_solve(varargin)
%{
    Compute Hx = b using LU factorization.

    Parameters
    --------------------
    Setting 1:
        H - hodlr
            Matrix in HODLR format - hodlr class.

        b - double/single


    Setting 2:
        l - hodlr
            Matrix in HODLR format - hodlr class.

        U - hodlr
            Matrix in HODLR format - hodlr class.
            
        b - double/single

    Returns
    --------------------
    x - double
        The solution.
%}

    if nargin == 2
        H = varargin{1}; 
        if ~(isa(H, 'hodlr') & ~isa(H, 'mphodlr'))
            error('Please ensure the first input is of a HODLR matrix.');
        end
        
        b = varargin{2}; 
        [L, U] = hlu(H, 'hodlr');
        y = htrsl(L, b);
        [m, n] = size(y);
        x = htrsu(U, y, 2);

    elseif nargin == 3
        L = varargin{1}; 
        U = varargin{2}; 

        if ~(isa(L, 'hodlr') & ~isa(L, 'mphodlr'))
            error('Please ensure the first input is of a HODLR matrix.');
        end

        if ~(isa(U, 'hodlr') & ~isa(U, 'mphodlr'))
            error('Please ensure the second input is of a HODLR matrix.');
        end
        
        b = varargin{3};
        y = htrsl(L, b);
        x = htrsu(U, y, 2);
    else
        error('Please enter correct number of inputs.');
    end
end