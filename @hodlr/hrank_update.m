function H = hrank_update(H, U, V, threshold)
    % Similar to hadd.m, but result in a HODLR format
    if isempty(H.D)
        [m, n, m1, m2, n1, n2] = hsize(H);
        [H.U1, H.V2] = compress_factors([H.U1, U(1:m1, :)], [H.V2', V(:, n1+1:end)']', threshold);
        [H.U2, H.V1] = compress_factors([H.U2, U(m1+1:end, :)], [H.V1', V(:, 1:n1)']', threshold);
        H.A11 = hrank_update(H.A11, U(1:m1,:), V(:,1:n1), threshold);
        H.A22 = hrank_update(H.A22, U(m1+1:end,:), V(:,n1+1:end), threshold);
    else
        H.D = H.D + U * V;
    end
end