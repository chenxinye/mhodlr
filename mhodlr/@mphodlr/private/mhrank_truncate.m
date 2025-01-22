function [Un, Vn] = mhrank_truncate(U, V, epsilon)
    [QU, RU] = qr(mchop(U), 0);
    [QV, RV] = qr(mchop(V'), 0);
    
    QU = mchop(QU);
    QV = mchop(QV);

    [U,S,V] = svd(full(mchop(RU * RV')));

    U = mchop(U);
    V = mchop(V);
    S = mchop(S);

    rnk = sum(abs(diag(S)) > S(1,1) * epsilon);
    Un = mchop(mchop(QU * U(:,1:rnk)) * S(1:rnk,1:rnk));
    Vn = mchop((QV * V(:,1:rnk))');
end