% A = load('data/bcspwr06.mat');
% A = A.Problem.A;
A = loadP64();

save('data/root_P64_cs128.mat', 'A')
%% mix precision
u1 = precision('s');
u2 = precision('h');
u3 = precision('b');
u4 = precision('q43');
u5 = precision('q52');

u_chain = prec_chain(u2, u3, u1, u4, u5);

epsilon = 1e-8;
aphA = amphodlr(u_chain, A, 5, 100, 'svd', epsilon); 
aprA = recover(aphA);
disp(aphA)
norm(aprA - A, 'fro') / norm(A, 'fro')