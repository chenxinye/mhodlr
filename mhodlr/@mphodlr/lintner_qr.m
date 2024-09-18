function [Q, R] = Lintner_QR(hA)
    B = hdot(hA.transpose(), hA);
    R =  hchol(B);
    Q = htrsu(hA, R, 1);
end
