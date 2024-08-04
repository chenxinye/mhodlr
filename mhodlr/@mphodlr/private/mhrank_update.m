function H = mhrank_update(H, U, V, epsilon)
    % Similar to hadd.m, but result in a HODLR format
    if isempty(H.D)
        [~, ~, m1, ~, n1, ~] = hsize(H);
        [H.U1, H.V2] = hrank_truncate([H.U1, U(1:m1, :)], [H.V2', V(:, n1+1:end)']', epsilon);
        [H.U2, H.V1] = hrank_truncate([H.U2, U(m1+1:end, :)], [H.V1', V(:, 1:n1)']', epsilon);
        H.A11 = mhrank_update(H.A11, U(1:m1,:), V(:,1:n1), epsilon);
        H.A22 = mhrank_update(H.A22, U(m1+1:end,:), V(:,n1+1:end), epsilon);
    else
        H.D = mchop(mchop(H.D) + mchop(U) * mchop(V));
    end
end