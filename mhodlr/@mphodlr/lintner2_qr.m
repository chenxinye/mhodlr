function [Q, R] = lintner2_qr(hA)
    B = hdot(hA.transpose(), hA);
    R = hchol(B);
    Q = htrsu(hA, R, 1);
    [Q, Rd] = lintner_qr(Q);
    R = hdot(Rd, R);
end
