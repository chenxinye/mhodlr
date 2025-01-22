
function [Q, R] = mhqr(H, method, prec)
%{
    The function is to perform QR factorization based on HODLR representation.

    Parameters
    --------------------
    A - hodlr
        Input matrix - hodlr class / dense tyle.

    method - str
        Options: 'lintner', 'bebendorf' and 'kressner'.
    
    prec - precision
        Precision to simulate the factorization.

    Returns
    --------------------
    Q, R - hodlr 
        Return matrix in hodlr class.
%}

    set_prec(prec);
    if strcmp(method, "lintner")
        [Q, R] = mlintner_qr(H);

    elseif strcmp(method, "bebendorf")
        [Q, R] = bebendorf_qr(H);

    elseif strcmp(method, "kressner")
        [Q, R] = mkressner_qr(H);
    end
end
    
    