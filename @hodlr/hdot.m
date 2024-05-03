function C = hdot(A, B, varargin) 
%{
    Compute dot product of A and B.

    Parameters
    --------------------
    A - hodlr | double
        Matrix in HODLR format or double array.
    
    B - hodlr | double
        Matrix in HODLR format or double array.

    oformat - str, default='hodlr'
    Output format
    


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
            oformat = 'double';
        end
    end

    C = hdot_double(A, B); 
    
    if strcmp(oformat, 'hodlr')
        C = hodlr(C);
    end
    
end