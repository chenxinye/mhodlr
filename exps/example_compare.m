
A = rand(5000, 5000);
x = rand(5000, 1);


tic;
b = A * x;
time = toc;

depth = 99; min_block_size = 200; epsilon = 1e-8;
method = 'svd';

hodlroption ('compression', method);
hodlroption('block-size', min_block_size);
hodlroption('threshold', epsilon);
addpath("mhodlr")


tic;
H = hodlr(A);
time_construct1 = toc;

tic;
b = H * x;
time1 = toc;

tic;
H = hodlr(A, depth, min_block_size, method, epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default
time_construct2 = toc;

tic;
b = hdot(H, x, 'dense');
time2 = toc;
