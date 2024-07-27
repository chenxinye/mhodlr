%% matrix 1
clear all

A = load('data/3-5000/root_P64_cs128.mat');
A =  A.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

depths = [2, 5, 8];
vareps = [1e-14, 1e-12, 1e-10, 1e-08, 1e-06, 1e-04, 1e-02];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_back_list1  = zeros(n_eps, n_d);
err_back_list2  = zeros(n_eps, n_d);
err_back_list3  = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        [L1, U1] = hlu(aphA, 'hodlr');
        [L2, U2] = mhlu(aphA, u2, 'hodlr');
        [L3, U3] = mhlu(aphA, u4, 'hodlr');

        recover_LU1 = hdot(L1, U1, 'dense');
        recover_LU2 = hdot(L2, U2, 'dense');
        recover_LU3 = hdot(L3, U3, 'dense');

        err_back1 = norm(A - recover_LU1, 'fro') / norm(A, 'fro');
        err_back2 = norm(A - recover_LU2, 'fro') / norm(A, 'fro');
        err_back3 = norm(A - recover_LU3, 'fro') / norm(A, 'fro');

        err_back_list1(i, j) = err_back1;
        err_back_list2(i, j) = err_back2;
        err_back_list3(i, j) = err_back3;
    end
end


save("results/lu1_depth2.mat", 'err_back_list1');
save("results/lu1_depth5.mat", 'err_back_list2');
save("results/lu1_depth8.mat", 'err_back_list3');



%% matrix 2
clear all

A = load('data/3-5000/cavity18.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

depths = [2, 5, 8];
vareps = [1e-14, 1e-12, 1e-10, 1e-08, 1e-06, 1e-04, 1e-02];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_back_list1  = zeros(n_eps, n_d);
err_back_list2  = zeros(n_eps, n_d);
err_back_list3  = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        [L1, U1] = hlu(aphA, 'hodlr');
        [L2, U2] = mhlu(aphA, u2, 'hodlr');
        [L3, U3] = mhlu(aphA, u4, 'hodlr');

        recover_LU1 = hdot(L1, U1, 'dense');
        recover_LU2 = hdot(L2, U2, 'dense');
        recover_LU3 = hdot(L3, U3, 'dense');

        err_back1 = norm(A - recover_LU1, 'fro') / norm(A, 'fro');
        err_back2 = norm(A - recover_LU2, 'fro') / norm(A, 'fro');
        err_back3 = norm(A - recover_LU3, 'fro') / norm(A, 'fro');

        err_back_list1(i, j) = err_back1;
        err_back_list2(i, j) = err_back2;
        err_back_list3(i, j) = err_back3;
    end
end




save("results/lu2_depth2.mat", 'err_back_list1');
save("results/lu2_depth5.mat", 'err_back_list2');
save("results/lu2_depth8.mat", 'err_back_list3');





%% matrix 3
clear all

A = load('data/3-5000/ex37.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

depths = [2, 5, 8];
vareps = [1e-14, 1e-12, 1e-10, 1e-08, 1e-06, 1e-04, 1e-02];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_back_list1  = zeros(n_eps, n_d);
err_back_list2  = zeros(n_eps, n_d);
err_back_list3  = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        [L1, U1] = hlu(aphA, 'hodlr');
        [L2, U2] = mhlu(aphA, u2, 'hodlr');
        [L3, U3] = mhlu(aphA, u4, 'hodlr');

        recover_LU1 = hdot(L1, U1, 'dense');
        recover_LU2 = hdot(L2, U2, 'dense');
        recover_LU3 = hdot(L3, U3, 'dense');

        err_back1 = norm(A - recover_LU1, 'fro') / norm(A, 'fro');
        err_back2 = norm(A - recover_LU2, 'fro') / norm(A, 'fro');
        err_back3 = norm(A - recover_LU3, 'fro') / norm(A, 'fro');

        err_back_list1(i, j) = err_back1;
        err_back_list2(i, j) = err_back2;
        err_back_list3(i, j) = err_back3;
    end
end



save("results/lu3_depth2.mat", 'err_back_list1');
save("results/lu3_depth5.mat", 'err_back_list2');
save("results/lu3_depth8.mat", 'err_back_list3');


%% matrix 4
clear all

A = load('data/3-5000/psmigr_1.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

depths = [2, 5, 8];
vareps = [1e-14, 1e-12, 1e-10, 1e-08, 1e-06, 1e-04, 1e-02];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_back_list1  = zeros(n_eps, n_d);
err_back_list2  = zeros(n_eps, n_d);
err_back_list3  = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        [L1, U1] = hlu(aphA, 'hodlr');
        [L2, U2] = mhlu(aphA, u2, 'hodlr');
        [L3, U3] = mhlu(aphA, u4, 'hodlr');

        recover_LU1 = hdot(L1, U1, 'dense');
        recover_LU2 = hdot(L2, U2, 'dense');
        recover_LU3 = hdot(L3, U3, 'dense');

        err_back1 = norm(A - recover_LU1, 'fro') / norm(A, 'fro');
        err_back2 = norm(A - recover_LU2, 'fro') / norm(A, 'fro');
        err_back3 = norm(A - recover_LU3, 'fro') / norm(A, 'fro');

        err_back_list1(i, j) = err_back1;
        err_back_list2(i, j) = err_back2;
        err_back_list3(i, j) = err_back3;
    end
end


save("results/lu4_depth2.mat", 'err_back_list1');
save("results/lu4_depth5.mat", 'err_back_list2');
save("results/lu4_depth8.mat", 'err_back_list3');