function [Q, R] = mlintner2_qr(hA)
    B = mhdot(hA.transpose(), hA);
    R = mhchol(B);
    Q = mhtrsu(hA, R, 1);
    [Q, Rd] = mlintner_qr(Q);
    Q = hmchop(Q)
    R = mhdot(Rd, R);
end
