function [Q, R] = hqr(H, method)
%{
    The function is used for the operation of summation or subtraction for HODLR matrix A and B.

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

    if strcmp(method, 'lintner')
        [Q, R] = lintner_qr(H);

    elseif strcmp(method, 'bebendorf')
        [Q, R] = bebendorf_qr(H);

    elseif strcmp(method, 'kressner')
        [Q, R] = kressner_qr(H);
    end
end

