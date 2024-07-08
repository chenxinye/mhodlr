A = load('data/lshp3025.mat');
A =  A.Problem.A;
%% mix precision
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u2, u3, u4, u5);

epsilon = 1e-2;
depth = 5;
aphA = amphodlr(u_chain, A, depth, 10, 'svd', epsilon); 
aprA = recover(aphA);


disp(aphA)

disp(epsilon)
norm(aprA - A, 'fro') / norm(A, 'fro')
bound_err = (2 * sqrt(2 * aphA.bottom_level) + 1) * epsilon 

VA = plot_hmat_prec(aphA, A, aphA.precIndex,1);

