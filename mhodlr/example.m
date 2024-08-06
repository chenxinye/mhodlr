rng(0);
A = rand(100);
x = rand(100, 1); % define vector
X = rand(100, 80); % define matrix

% Usual call for full working precision 
hA = hodlr(A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
rA = recover(hA);
norm(rA - A, 2)  % Compute the error


b = hdot(hA, x);
err = norm(recover(b) - A * x, 'fro');
disp(err);


b = hdot(hA, x, 'dense');
err = norm(b - A * x, 'fro');
disp(err);

B = hdot(hA, X);
err = norm(recover(B) - A * X, 'fro');
disp(err);

B = hdot(hA, X, 'dense');
err = norm(B - A * X, 'fro');
disp(err);



[L, U] = hlu(hA);
err = norm(hdot(L, U, 'dense') - A, 'fro');
disp(err);

[L, U] = hlu(hA, 'dense');
err = norm(L * U - A, 'fro');
disp(err);

