A = load('data/bcspwr09.mat');
A = A.Problem.A;
%% mix precision
u1 = precision('s');
u2 = precision('h');
u3 = precision('q52');

u_chain = prec_chain(u1, u3, u2);

epsilon = 1e-5;
depth = 5;
aphA = amphodlr(u_chain, A, depth, 10, 'svd', epsilon); 
aprA = recover(aphA);


disp(aphA)

disp(epsilon)
norm(aprA - A, 'fro') / norm(A, 'fro')
bound_err = (2 * sqrt(2 * aphA.bottom_level) + 1) * epsilon 