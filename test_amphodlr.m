
%% mix precision
u1 = precision('q43');
u2 = precision('q52');
u3 = precision('b');
u4 = precision('h');
u5 = precision('s');
u_chain = prec_chain(u1, u2, u3, u4, u5);
disp('Precision [q43, q52, b, h, s], eps=0.001')
epsilon = 0.001;

n = 100;
A = rand(n, n);
aphA = amphodlr(u_chain, A, 5, 2, 'svd', epsilon);
aprA = recover(aphA);
norm(aprA - A, 'fro')

phA = mphodlr(u_chain, A, 5, 20, 'svd', epsilon);
prA = recover(phA);
norm(prA - A, 'fro')

A = spdiags(ones(n, 1) * [2 3 -1 0.5],  -1:1, n, n); 

aphA = amphodlr(u_chain, A, 5, 2, 'svd', epsilon);
aprA = recover(aphA);
norm(aprA - A, 'fro')

phA = mphodlr(u_chain, A, 5, 20, 'svd', epsilon);
prA = recover(phA);
norm(prA - A, 'fro')
