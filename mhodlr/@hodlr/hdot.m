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

    
    if strcmp(oformat, 'hodlr') 
        C = hdot_hodlr(A, B); 
        % C = hodlr(C, A.max_level, A.min_block_size, A.method);
    else
        C = hdot_dense(A, B); 
    end

end