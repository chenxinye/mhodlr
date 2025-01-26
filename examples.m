
addpath("./mhodlr/")
rng(0); %fix randomness
A = rand(20, 20);
depth = 99;
min_block_size = 2;
epsilon = 1e-14;
hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default

disp([norm(A, 'fro'), hnorm(hA, 'fro')])

hA = htruncate(hA, 10);
rA = recover(hA);
disp(norm(rA - A, 2));

disp(hA);
disp(hA.A11);

%%% Matrix transpose and inverse

rng(1); %fix randomness
A = rand(50);
depth = 5;
min_block_size = 2;
epsilon = 1e-12;
hA = hodlr(A, depth, min_block_size, 'svd', epsilon, 20); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default
iA = inverse(hA); 
disp(norm(hdot(hA, iA, 'dense') - eye(50), 'fro'))

iA = inverse(hA, 'dense');  % return an the inverse in dense format
disp(norm(iA * A - eye(50), 'fro'))

disp(norm(hadd(hA.transpose(), A', '-', 'dense'), 'fro'))

% mixed precision
u = precision('s');
set_prec(u);
iA = minverse(hA); 
disp(norm(hdot(hA, iA, 'dense') - eye(50), 'fro'))

iA = minverse(hA, 'dense');  % return an the inverse in dense format
disp(norm(iA * A - eye(50), 'fro'))



%%% Matrix summation and subtraction
rng(0); %fix randomness
A = rand(500);
B = rand(500);

depth = 5;
min_block_size = 2;
epsilon = 1e-8;
hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default
hB = hodlr(B, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default

disp(norm(hadd(hA, B, '-', 'dense') - (A-B), 'fro'))
disp(norm(hadd(hA, hB, '-', 'dense') - (A-B), 'fro'))
disp(norm(hadd(A, hB, '-', 'dense') - (A-B), 'fro'))
disp(norm(recover(hadd(hA, B, '-', 'hodlr')) - (A-B), 'fro'))
disp(norm(recover(hadd(hA, hB, '-', 'hodlr')) - (A-B), 'fro'))
disp(norm(recover(hadd(A, hB, '-', 'hodlr')) - (A-B), 'fro'))

disp(norm(hadd(hA, B, '+', 'dense') - (A+B), 'fro'))
disp(norm(hadd(hA, hB, '+', 'dense') - (A+B), 'fro'))
disp(norm(hadd(A, hB, '+', 'dense') - (A+B), 'fro'))
disp(norm(recover(hadd(hA, B, '+', 'hodlr')) - (A+B), 'fro'))
disp(norm(recover(hadd(hA, hB, '+', 'hodlr')) - (A+B), 'fro'))
disp(norm(recover(hadd(A, hB, '+', 'hodlr')) - (A+B), 'fro'))

u = precision('h');
set_prec(u);
C1 = add(hA, hB, hA);
C2 = A + B + A;
disp(norm(C1.dense - C2, 'fro'));

C1 = sub(hA, hB, hA);
C2 = A - B - A;
disp(norm(C1.dense - C2, 'fro'));

C1 = madd(hA, hB, hA);
C2 = A + B + A;
disp(norm(C1.dense - C2, 'fro'));

C1 = msub(hA, hB, hA);
C2 = A - B - A;
disp(norm(C1.dense - C2, 'fro'));


C1 = madd(hA, B, hA);
C2 = A + B + A;
disp(norm(C1.dense - C2, 'fro'));

C1 = msub(hA, B, hA);
C2 = A - B - A;
disp(norm(C1.dense - C2, 'fro'));

%%% Build/Recover/Product
rng(0);
A = rand(100);
x = rand(100, 2); % define vector
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

B = dot(hA, hA, hA);
err = norm(B.dense - A * A * A, 'fro');
disp(err);



%%% H-LU factorization
[L, U] = hlu(hA);
err = norm(hdot(L, U, 'dense') - A, 'fro');
disp(err);

[L, U] = hlu(hA, 'dense');
err = norm(L * U - A, 'fro');
disp(err);

[L, U] = hlu(hA, 'dense');
err = norm(L * U - A, 'fro');
disp(err);


u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');

u_chain = prec_chain(u1, u2, u3, u4);
depth=5;
aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
mphA = mphodlr(u_chain, A, depth, 10, 'svd', eps); 

u = precision('h');
set_prec(u);
[L, U] = mhlu(hA, 'hodlr');
err = norm(hdot(L, U, 'dense') - A, 'fro');
disp(err);

u = precision('b');
set_prec(u);
[L, U] = mhlu(mphA, 'hodlr');
err = norm(hdot(L, U, 'dense') - A, 'fro');
disp(err);

u = precision('s');
set_prec(u);
[L, U] = mhlu(aphA, 'hodlr');
err = norm(hdot(L, U, 'dense') - A, 'fro');
disp(err);




%%% H-Cholesky factorization

rng(0);
R = rand(80);
A = R'*R; % Generate symmetric positive definite matrix

% Usual call for full working precision 
hA = hodlr(A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.

R = hchol(hA); % return a block upper triangular HODLR matrix

disp(norm(hdot(R.transpose(), R, 'dense') - A, 'fro'))

R = hchol(hA, 'dense'); % return a upper triangulr matrix
norm(R'*R - A, 'fro')

% Create precisions for each level; Level 1 use precision u1, level 2 use precision u2, ...
u1 = precision('q43');
u2 = precision('q52');
u3 = precision('b');
u4 = precision('s');
u_chain = prec_chain(u1, u2, u3, u4);


% Call mixed precision HODLR representation
amphA = amphodlr(u_chain, A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
amprA = recover(amphA);
norm(amprA - A,2) % Compute the error

set_prec(u4);
R = mhchol(amphA); % or R = mhchol(hA, u4);
disp(norm(hdot(R.transpose(), R, 'dense') - A, 'fro'))



%%% QR factorization
rng(0);
A = rand(100);
depth = 5;
min_block_size = 2;
epsilon = 1e-8;
hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default

[Q, R] = hqr(hA, 'lintner');

disp(norm(hdot(Q.transpose(), Q, 'dense') - eye(100), 'fro'))
disp(norm(hdot(Q, R, 'dense') - A, 'fro'))



%%% Generate special matrices
u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');

u_chain = prec_chain(u1, u2, u3, u4);
A = hodlr('eye', 100, 10);
B = mphodlr('eye', 100, 10);
C = amphodlr('eye', 100, 10);

disp(norm(recover(A) - eye(100), 'fro'))
disp(norm(recover(B) - eye(100), 'fro'))
disp(norm(recover(C) - eye(100), 'fro'))

A = hodlr('ones', 100, 10);
B = mphodlr('ones', 100, 10);
C = amphodlr('ones', 100, 10);

disp(norm(recover(A) - ones(100), 'fro'))
disp(norm(recover(B) - ones(100), 'fro'))
disp(norm(recover(C) - ones(100), 'fro'))



% QR 
rng(0); %fix randomness
A = rand(60, 50);
depth = 88;
min_block_size = 2;
epsilon = 1e-14;
hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default

[Q, R] = hqr(hA, 'dk');

disp(norm(hdot(Q,R).dense -A))

rng(0); %fix randomness
A = rand(30, 30);
depth = 5;
min_block_size = 2;
epsilon = 1e-14;
hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default

[Q, R] = lintner_qr(hA);
norm(Q.dense*R.dense - A)
norm(Q.transpose.dense*Q.dense - eye(30))

[Q, R] = lintner2_qr(hA);
norm(Q.dense*R.dense - A)
norm(Q.transpose.dense*Q.dense - eye(30))


% multiple precision QR 

u = precision('s');
set_prec(u);
[Q, R] = mhqr(hA, 'lt');
[Q, R] = mhqr(hA, 'lt2');
