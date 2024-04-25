n = 10;
A = spdiags(ones(n, 1) * [1 3  -1],  -1:1, n, n);
B = spdiags(ones(n, 1) * [1 3  1],  -1:1, n, n);

[U, V] = compress_test(A(1:5, 5:end), ...
    'svd', 1e-8);
norm(U * V - A(1:5, 5:end), 2);

hA = hodlr(A);

hB = hodlr(B);
R = hchol(hB);
norm(full(R' * R - B),2)
recover(full(hB))
full(R' * R)

hC1 = hadd(hA, A, '+')
hC2 = hadd(hA, hA, '+')
hC3 = hadd(hA, hA, '+')
hC4 = hadd(A, B, '+')

thred = 10^(-10);
% [Un,Vn] = compress_factors(U, V, thred)


[L, U] = hlu(hA, 10^(-10));

norm(full(L * U - A),2)

invhA = hA.inverse_double()
invhA = hA.inverse_hodlr()
norm(full(hdot(invhA, A, 'double') - eye(n)),2)

[m,n] = hsize(hA)
[m] = hsize(hA)
[m, n, a, b, c, d] = hsize(hA)

C = hdot(hA, A);
norm(C - A*A,2)


C = hdot(A, hA, 'double');
norm(C - A*A,2)

C = hdot(A, A, 'double');
norm(C - A*A,2)


C = recover(hdot(hA, hA));
norm(C - A*A,2)


bAT = hA.transpose();
norm(A.' - recover(bAT),2)

issqaure(hA);

% misc

% Interpolative decomposition
A = [34, 58, 52hA; 59, 89, 80; 17, 29, 26];
B = [58, 34; 89, 59; 29, 17];
C = [0, 1, 29/33; 1, 0, 1/33];