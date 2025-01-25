function [Q, R] = hqr(H, method)
%{
    The function is to perform QR factorization based on HODLR representation.

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

    if strcmp(method, 'lt')
        [Q, R] = lintner_qr(H);

    elseif strcmp(method, 'lt2')
        [Q, R] = lintner2_qr(H);

    elseif strcmp(method, 'dk')
        [Q, R] = kressner_qr(H);
    end
end

