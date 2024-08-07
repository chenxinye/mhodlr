function C = mhdot(A, B, prec, varargin) 
%{
    Compute dot product of A and B.

    Parameters
    --------------------
    A - hodlr | double
        Matrix in HODLR format or double array.
    
    B - hodlr | double
        Matrix in HODLR format or double array.

    prec - precision
        Precision for the matrix-vector product.

    oformat - str, default='hodlr'
        Output format: 'hodlr' or 'dense'.
    

    Returns
    --------------------
    C - hodlr | double
        The matrix of product. 
%}
 
    if nargin == 3 
        oformat = 'hodlr';
    else
        if strcmp(varargin{1}, 'hodlr')
            oformat = 'hodlr';
        else
            oformat = 'dense';
        end
    end
    
    set_prec(prec);
    
    C = mhdot_dense(A, B); 
    
    if strcmp(oformat, 'hodlr')
        C = hodlr(C);
    end
    
end