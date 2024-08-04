function [Q, R] = hqr(H)
    HH = hdot(H.transpose(), H, 'hodlr');
    R = hchol(HH);
    Q = htrsu(H, R, 0);
end
