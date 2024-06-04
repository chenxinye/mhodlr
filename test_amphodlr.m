
%% mix precision
u1 = precision('s');
u2 = precision('h');
u3 = precision('b');
u4 = precision('q43');
u5 = precision('q52');

u_chain = prec_chain(u2, u3, u1, u4, u5);

u_chain_2 = prec_chain(u2, u3, u1, u4, u5); % mix order

u_chain_inv = prec_chain(u5, u4, u3, u2, u1);
disp('Precision [s, h, b, q52, q43], eps=0.001')
epsilon = 1e-12;

n = 1000;
A = rand(n, n);

aphA = amphodlr(u_chain, A, 5, 20, 'svd', epsilon); 
aprA = recover(aphA);
norm(aprA - A, 'fro') / norm(A, 'fro')

phA = mphodlr(u_chain, A, 5, 20, 'svd', epsilon);
prA = recover(phA);
norm(prA - A, 'fro')

phA = mphodlr(u_chain_inv, A, 5, 20, 'svd', epsilon);
prA = recover(phA);
norm(prA - A, 'fro')


A = spdiags(ones(n, 1) * [1 -1 1],  -1:1, n, n); 

aphA = amphodlr(u_chain, A, 5, 20, 'svd', epsilon);
aprA = recover(aphA);
norm(aprA - A, 'fro')


phA = mphodlr(u_chain, A, 5, 20, 'svd', epsilon);
prA = recover(phA);
norm(prA - A, 'fro')


phA = mphodlr(u_chain_inv, A, 5, 20, 'svd', epsilon);
prA = recover(phA);
norm(prA - A, 'fro')
