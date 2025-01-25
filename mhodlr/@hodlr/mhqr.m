
function [Q, R] = mhqr(H, method)
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

    global opt;
    set_prec(opt);

    if strcmp(method, "lintner")
        [Q, R] = mlintner_qr(H);

    elseif strcmp(method, "bebendorf")
        [Q, R] = mbf_qr(H);

    elseif strcmp(method, "kressner")
        [Q, R] = mkressner_qr(H);
    end
end
    
    
