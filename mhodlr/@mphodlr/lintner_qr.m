function [Q, R] = lintner_qr(hA)
    B = hdot(hA.transpose(), hA);
    R =  hchol(B);
    Q = htrsu(hA, R, 1);
end
