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

err_bound_list1 = zeros(n_eps, n_d);
err_bound_list2 = zeros(n_eps, n_d);
err_bound_list3 = zeros(n_eps, n_d);

norm_A = norm(A, 'fro');
n_A = size(A, 1);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        [L1, U1] = hlu(aphA, 'hodlr');
        norm_L = norm(recover(L1), 'fro');
        norm_U = norm(recover(U1), 'fro');
        eb1 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        [L2, U2] = mhlu(aphA, u2, 'hodlr');
        norm_L = norm(recover(L2), 'fro');
        norm_U = norm(recover(U2), 'fro');
        eb2 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        [L3, U3] = mhlu(aphA, u4, 'hodlr');
        norm_L = norm(recover(L3), 'fro');
        norm_U = norm(recover(U3), 'fro');
        eb3 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        recover_LU1 = hdot(L1, U1, 'dense');
        recover_LU2 = hdot(L2, U2, 'dense');
        recover_LU3 = hdot(L3, U3, 'dense');
        
        err_back1 = norm(A - recover_LU1, 'fro') / norm_A;
        err_back2 = norm(A - recover_LU2, 'fro') / norm_A;
        err_back3 = norm(A - recover_LU3, 'fro') / norm_A;

        err_back_list1(i, j) = err_back1;
        err_back_list2(i, j) = err_back2;
        err_back_list3(i, j) = err_back3;

        err_bound_list1(i, j) = eb1;
        err_bound_list2(i, j) = eb2;
        err_bound_list3(i, j) = eb3;
    end
end

save("results/lu1_1.mat", 'err_back_list1');
save("results/lu1_2.mat", 'err_back_list2');
save("results/lu1_3.mat", 'err_back_list3');

save("results/lu1_1_bound.mat", 'err_bound_list1');
save("results/lu1_2_bound.mat", 'err_bound_list2');
save("results/lu1_3_bound.mat", 'err_bound_list3');



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

err_bound_list1 = zeros(n_eps, n_d);
err_bound_list2 = zeros(n_eps, n_d);
err_bound_list3 = zeros(n_eps, n_d);

norm_A = norm(A, 'fro');
n_A = size(A, 1);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        [L1, U1] = hlu(aphA, 'hodlr');
        norm_L = norm(recover(L1), 'fro');
        norm_U = norm(recover(U1), 'fro');
        eb1 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        [L2, U2] = mhlu(aphA, u2, 'hodlr');
        norm_L = norm(recover(L2), 'fro');
        norm_U = norm(recover(U2), 'fro');
        eb2 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        [L3, U3] = mhlu(aphA, u4, 'hodlr');
        norm_L = norm(recover(L3), 'fro');
        norm_U = norm(recover(U3), 'fro');
        eb3 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        recover_LU1 = hdot(L1, U1, 'dense');
        recover_LU2 = hdot(L2, U2, 'dense');
        recover_LU3 = hdot(L3, U3, 'dense');
        
        err_back1 = norm(A - recover_LU1, 'fro') / norm_A;
        err_back2 = norm(A - recover_LU2, 'fro') / norm_A;
        err_back3 = norm(A - recover_LU3, 'fro') / norm_A;

        err_back_list1(i, j) = err_back1;
        err_back_list2(i, j) = err_back2;
        err_back_list3(i, j) = err_back3;

        err_bound_list1(i, j) = eb1;
        err_bound_list2(i, j) = eb2;
        err_bound_list3(i, j) = eb3;
    end
end


save("results/lu2_depth2.mat", 'err_back_list1');
save("results/lu2_depth5.mat", 'err_back_list2');
save("results/lu2_depth8.mat", 'err_back_list3');

save("results/lu2_1_bound.mat", 'err_bound_list1');
save("results/lu2_2_bound.mat", 'err_bound_list2');
save("results/lu2_3_bound.mat", 'err_bound_list3');



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

err_bound_list1 = zeros(n_eps, n_d);
err_bound_list2 = zeros(n_eps, n_d);
err_bound_list3 = zeros(n_eps, n_d);

norm_A = norm(A, 'fro');
n_A = size(A, 1);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        [L1, U1] = hlu(aphA, 'hodlr');
        norm_L = norm(recover(L1), 'fro');
        norm_U = norm(recover(U1), 'fro');
        eb1 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        [L2, U2] = mhlu(aphA, u2, 'hodlr');
        norm_L = norm(recover(L2), 'fro');
        norm_U = norm(recover(U2), 'fro');
        eb2 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        [L3, U3] = mhlu(aphA, u4, 'hodlr');
        norm_L = norm(recover(L3), 'fro');
        norm_U = norm(recover(U3), 'fro');
        eb3 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        recover_LU1 = hdot(L1, U1, 'dense');
        recover_LU2 = hdot(L2, U2, 'dense');
        recover_LU3 = hdot(L3, U3, 'dense');
        
        err_back1 = norm(A - recover_LU1, 'fro') / norm_A;
        err_back2 = norm(A - recover_LU2, 'fro') / norm_A;
        err_back3 = norm(A - recover_LU3, 'fro') / norm_A;

        err_back_list1(i, j) = err_back1;
        err_back_list2(i, j) = err_back2;
        err_back_list3(i, j) = err_back3;

        err_bound_list1(i, j) = eb1;
        err_bound_list2(i, j) = eb2;
        err_bound_list3(i, j) = eb3;
    end
end


save("results/lu3_depth2.mat", 'err_back_list1');
save("results/lu3_depth5.mat", 'err_back_list2');
save("results/lu3_depth8.mat", 'err_back_list3');

save("results/lu3_1_bound.mat", 'err_bound_list1');
save("results/lu3_2_bound.mat", 'err_bound_list2');
save("results/lu3_3_bound.mat", 'err_bound_list3');

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

err_bound_list1 = zeros(n_eps, n_d);
err_bound_list2 = zeros(n_eps, n_d);
err_bound_list3 = zeros(n_eps, n_d);

norm_A = norm(A, 'fro');
n_A = size(A, 1);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        [L1, U1] = hlu(aphA, 'hodlr');
        norm_L = norm(recover(L1), 'fro');
        norm_U = norm(recover(U1), 'fro');
        eb1 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        [L2, U2] = mhlu(aphA, u2, 'hodlr');
        norm_L = norm(recover(L2), 'fro');
        norm_U = norm(recover(U2), 'fro');
        eb2 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        [L3, U3] = mhlu(aphA, u4, 'hodlr');
        norm_L = norm(recover(L3), 'fro');
        norm_U = norm(recover(U3), 'fro');
        eb3 = lu_err_bound2(depth, eps, norm_A, norm_L, norm_U)/ norm_A;

        recover_LU1 = hdot(L1, U1, 'dense');
        recover_LU2 = hdot(L2, U2, 'dense');
        recover_LU3 = hdot(L3, U3, 'dense');
        
        err_back1 = norm(A - recover_LU1, 'fro') / norm_A;
        err_back2 = norm(A - recover_LU2, 'fro') / norm_A;
        err_back3 = norm(A - recover_LU3, 'fro') / norm_A;

        err_back_list1(i, j) = err_back1;
        err_back_list2(i, j) = err_back2;
        err_back_list3(i, j) = err_back3;

        err_bound_list1(i, j) = eb1;
        err_bound_list2(i, j) = eb2;
        err_bound_list3(i, j) = eb3;
    end
end


save("results/lu4_depth2.mat", 'err_back_list1');
save("results/lu4_depth5.mat", 'err_back_list2');
save("results/lu4_depth8.mat", 'err_back_list3');

save("results/lu4_1_bound.mat", 'err_bound_list1');
save("results/lu4_2_bound.mat", 'err_bound_list2');
save("results/lu4_3_bound.mat", 'err_bound_list3');
