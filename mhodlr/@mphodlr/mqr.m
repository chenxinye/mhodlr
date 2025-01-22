
function [Q, R] = hqr(H, method, prec)
%{
    The function is to perform QR factorization based on HODLR representation.

    Parameters
    --------------------
    A - hodlr
        Input matrix - hodlr class / dense tyle.

    method - str
        Options: 'lintner', 'bebendorf' and 'kressner'.
    
    Returns
    --------------------
    Q, R - hodlr 
        Return matrix in hodlr class.
%}

    set_prec(prec);
    if strcmp(method, 'lintner')
        [Q, R] = lintner_qr(H);

    elseif strcmp(method, 'bebendorf')
        [Q, R] = bebendorf_qr(H);

    elseif strcmp(method, 'kressner')
        [Q, R] = mkressner_qr(H);
    end
end
    
    