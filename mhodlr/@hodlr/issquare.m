function ch = issquare(H)
%{
    The function is to check whether HODLR matrix H is square or not.

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
    

    Returns
    --------------------
    ch - boolean
        `1` indicates H is square matrix, `0` indicates otherwise. 
%}
 
 
    if strcmp(class(H), 'hodlr') | strcmp(class(H), 'mphodlr') | strcmp(class(H), 'amphodlr')
        if isempty(H.D)
            su1 = size(H.U1, 1);
            su2 = size(H.U2, 1);
            sv1 = size(H.V1, 2);
            sv2 = size(H.V2, 2);
            
            m = su1 + su2;
            n = sv1 + sv2;
            ch = m == n;
        else
            [m, n] = size(H.D);
            ch = m == n;
        end

            
    elseif class(H) == 'double' or 'single'
        [m, n] = size(H);
        ch = m == n;
    else
        error('Please enter correct type of input.');
    end
end