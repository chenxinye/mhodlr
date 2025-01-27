function [Q, R] = mhqr(H, method)
%{
    The function is to perform QR factorization in multiple precision based on HODLR representation.

    Parameters
    --------------------
    A - hodlr
        Input matrix - hodlr class / dense tyle.

    method - str
        Options: 'lt', 'lt2' and 'dk'.
    
    Returns
    --------------------
    Q, R - hodlr 
        Return matrix in hodlr class.
%}

    global opt;
    set_prec(opt);
    
    if strcmp(method, 'lt')
        [Q, R] = mlintner_qr(H);
    
    elseif strcmp(method, 'lt2')
        [Q, R] = mlintner2_qr(H);

    elseif strcmp(method, 'kd')
        [Q, R] = mkressner_qr(H);
    end
end
    
