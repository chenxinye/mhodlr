% P64
clear all
A = load('data/3-5000/root_P64_cs128.mat');
A =  A.A;
disp(size(A));
u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

vareps = [10e-9, 10e-08, 10e-07, 10e-06, 10e-5, 10e-04, 10e-03, 10e-02]; % u2.u, u3.u
depths = [2,  5,  8];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_list = zeros(n_eps, n_d);
err_list_amp = zeros(n_eps, n_d);
err_bound_list = zeros(n_eps, n_d);

for j = 1:n_d
    d = depths(j);
    for i = 1:n_eps
        eps = vareps(i);

        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps); 
        aprA = recover(aphA);
    
        err_list(i, j) = norm(rA - A, 'fro') / norm(A, 'fro');
        err_list_amp(i, j) = norm(aprA - A, 'fro') / norm(A, 'fro');
        err_bound_list(i, j) = err_bound(aphA.bottom_level, eps);
    
        disp([err_list(i, j), err_list_amp(i, j), err_bound_list(i, j)]);
        disp(aphA.precIndex); 
        disp('--------------------------'); 
    end 
    disp('*******************************'); 
end

save("results/rcerr_err_list1.mat", 'err_list');
save("results/rcerr_err_list_amp1.mat", 'err_list_amp');
save("results/rcerr_err_bound_list1.mat", 'err_bound_list');



% ex37
clear all
A = load('data/3-5000/ex37.mat');
A =  A.Problem.A;
disp(size(A));
u1 = precision('d')
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

vareps = [10e-9, 10e-08, 10e-07, 10e-06, 10e-5, 10e-04, 10e-03, 10e-02]; % u2.u, u3.u
depths = [2,  5,  8];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_list = zeros(n_eps, n_d);
err_list_amp = zeros(n_eps, n_d);
err_bound_list = zeros(n_eps, n_d);

for j = 1:n_d
    d = depths(j);
    for i = 1:n_eps
        eps = vareps(i);

        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps); 
        aprA = recover(aphA);
    
        err_list(i, j) = norm(rA - A, 'fro') / norm(A, 'fro');
        err_list_amp(i, j) = norm(aprA - A, 'fro') / norm(A, 'fro');
        err_bound_list(i, j) = err_bound(aphA.bottom_level, eps);
    
        disp([err_list(i, j), err_list_amp(i, j), err_bound_list(i, j)]);
        disp(aphA.precIndex); 
        disp('--------------------------'); 
    end 
    disp('*******************************'); 
end


save("results/rcerr_err_list2.mat", 'err_list');
save("results/rcerr_err_list_amp2.mat", 'err_list_amp');
save("results/rcerr_err_bound_list2.mat", 'err_bound_list');



% 1138_bus
clear all
A = load('data/1138_bus.mat');
A =  A.Problem.A;
disp(size(A));
u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

vareps = [10e-9, 10e-08, 10e-07, 10e-06, 10e-5, 10e-04, 10e-03, 10e-02]; % u2.u, u3.u
depths = [2,  5,  8];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_list = zeros(n_eps, n_d);
err_list_amp = zeros(n_eps, n_d);
err_bound_list = zeros(n_eps, n_d);

for j = 1:n_d
    d = depths(j);
    for i = 1:n_eps
        eps = vareps(i);

        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps); 
        aprA = recover(aphA);
    
        err_list(i, j) = norm(rA - A, 'fro') / norm(A, 'fro');
        err_list_amp(i, j) = norm(aprA - A, 'fro') / norm(A, 'fro');
        err_bound_list(i, j) = err_bound(aphA.bottom_level, eps);
    
        disp([err_list(i, j), err_list_amp(i, j), err_bound_list(i, j)]);
        disp(aphA.precIndex); 
        disp('--------------------------'); 
    end 
    disp('*******************************'); 
end


save("results/rcerr_err_list3.mat", 'err_list');
save("results/rcerr_err_list_amp3.mat", 'err_list_amp');
save("results/rcerr_err_bound_list3.mat", 'err_bound_list');

