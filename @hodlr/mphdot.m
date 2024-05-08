function C = mphdot(A, B, varargin) 
%{
    Compute dot product of A and B in mixed precision.

    Parameters
    --------------------
    A - hodlr | double
        Matrix in HODLR format or dense array.
    
    B - hodlr | double
        Matrix in HODLR format or dense array.

    oformat - str, default='hodlr'
        Output format: 'hodlr' or 'dense'.

    Returns
    --------------------
    C - hodlr | double
        The matrix of product. 
%}
 
    if nargin == 2 
        oformat = 'hodlr';
    else
        if strcmp(varargin{1}, 'hodlr')
            oformat = 'hodlr';
        else
            oformat = 'dense';
        end
    end

    C = mphdot_dense(A, B); 
    
    if strcmp(oformat, 'hodlr')
        C = hodlr(C);
    end
    
end