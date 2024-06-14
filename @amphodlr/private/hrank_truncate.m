function [Un, Vn] = hrank_truncate(U, V, epsilon)
    [QU, RU] = qr(U, 0);
    [QV, RV] = qr(V', 0);
    
    [U,S,V] = svd(full(RU * RV'));

    i = min(size(S));
    rnk = i-length(find(diag(S(1:i,1:i)) <= S(1,1)*epsilon));
    Un = QU * U(:,1:rnk) * S(1:rnk,1:rnk);
    Vn = (QV * V(:,1:rnk))';
end