function [Un, Vn] = compress_factors(U, V, epsilon)
    [QU, RU] = qr(U, 0);
    [QV, RV] = qr(V', 0);
    
    [U,S,V] = svd(full(RU * RV'));

    rnk = sum(abs(diag(S)) > S(1,1) * epsilon);
    Un = QU * U(:,1:rnk) * S(1:rnk,1:rnk);
    Vn = (QV * V(:,1:rnk))';
end