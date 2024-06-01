function H = hrank_update(H, U, V, epsilon)
    % Similar to hadd.m, but result in a HODLR format
    if isempty(H.D)
        [m, n, m1, m2, n1, n2] = hsize(H);
        [H.U1, H.V2] = compress_factors([H.U1, U(1:m1, :)], [H.V2', V(:, n1+1:end)']', epsilon);
        [H.U2, H.V1] = compress_factors([H.U2, U(m1+1:end, :)], [H.V1', V(:, 1:n1)']', epsilon);
        H.A11 = hrank_update(H.A11, U(1:m1,:), V(:,1:n1), epsilon);
        H.A22 = hrank_update(H.A22, U(m1+1:end,:), V(:,n1+1:end), epsilon);
    else
        H.D = H.D + U * V;
    end
end