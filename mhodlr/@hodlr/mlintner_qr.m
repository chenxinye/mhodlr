function [Q, R] = mlintner_qr(hA)
    global opt;
    B = hmchop(hdot(hA.transpose(), hA));
    R = mhchol(B, opt);
    Q = mhtrsu(hA, R, 1);
end
