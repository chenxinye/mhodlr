%% kernel matrix 1
rng(0)

n_sample = 10;
x = rand(1, 2000);
y = rand(1, 2000);
kernel_mat = kernel1(x, y);

rng(42);
v =  unifrnd(-1, 1, n_sample, 2000);

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
ref_err_back_list  = zeros(n_eps, n_d);
err_bound_list = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, kernel_mat, depth, 10, 'svd', eps); 

        for k=1:n_sample
            x = v(k, :)';
            
            hb1 = mhdot(aphA, x, u2, 'dense');
            hb2 = mhdot(aphA, x, u4, 'dense');
            rhb = hdot(aphA, x, 'dense');
            
            b = kernel_mat * x;
            
            err_back1 = norm(b - hb1, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            err_back2 = norm(b - hb2, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            ref_err_back = norm(b - rhb, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            err_bound = 2*(sqrt(2)+1)*sqrt(2^(depth + 1)+2^(depth-1))*eps;

            err_back_list1(i, j) = err_back_list1(i, j) + err_back1;
            err_back_list2(i, j) = err_back_list2(i, j) + err_back2;
            ref_err_back_list(i, j) = ref_err_back_list(i, j) + ref_err_back;
            err_bound_list(i, j) = err_bound_list(i, j) + err_bound;
        end
    end
end

save("results/prod1_err_back1.mat", 'err_back_list1');
save("results/prod1_err_back2.mat", 'err_back_list2');
save("results/prod1_ref_err_back.mat", 'ref_err_back_list');
save("results/prod1_bound.mat", 'err_bound_list');


%% kernel matrix 2
rng(0)

n_sample = 10;
x = unifrnd(-1, 1, 2000, 2);
y = unifrnd(-1, 1, 2000, 2);
kernel_mat = kernel2(x, y);

rng(42);
v =  unifrnd(-1, 1, n_sample, 2000);

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
ref_err_back_list  = zeros(n_eps, n_d);
err_bound_list = zeros(n_eps, n_d);


for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, kernel_mat, depth, 10, 'svd', eps); 

        for k=1:n_sample
            x = v(k, :)';
            
            hb1 = mhdot(aphA, x, u2, 'dense');
            hb2 = mhdot(aphA, x, u4, 'dense');
            rhb = hdot(aphA, x, 'dense');
            
            b = kernel_mat * x;
            
            err_back1 = norm(b - hb1, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            err_back2 = norm(b - hb2, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            ref_err_back = norm(b - rhb, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            err_bound = 2*(sqrt(2)+1)*sqrt(2^(depth + 1)+2^(depth-1))*eps;

            err_back_list1(i, j) = err_back_list1(i, j) + err_back1;
            err_back_list2(i, j) = err_back_list2(i, j) + err_back2;
            ref_err_back_list(i, j) = ref_err_back_list(i, j) + ref_err_back;
            err_bound_list(i, j) = err_bound_list(i, j) + err_bound;
        end
    end
end

save("results/prod2_err_back1.mat", 'err_back_list1');
save("results/prod2_err_back2.mat", 'err_back_list2');
save("results/prod2_ref_err_back.mat", 'ref_err_back_list');
save("results/prod2_bound.mat", 'err_bound_list');


%% kernel matrix 3
clear all
rng(0)

n_sample = 10;
x = unifrnd(-1, 1, 2000, 2);
y = unifrnd(-1, 1, 2000, 2);
kernel_mat = kernel3(x, y);

rng(42);
v =  unifrnd(-1, 1, n_sample, 2000);

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
ref_err_back_list  = zeros(n_eps, n_d);
err_bound_list = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, kernel_mat, depth, 10, 'svd', eps); 

        for k=1:n_sample
            x = v(k, :)';
            
            hb1 = mhdot(aphA, x, u2, 'dense');
            hb2 = mhdot(aphA, x, u4, 'dense');
            rhb = hdot(aphA, x, 'dense');
            
            b = kernel_mat * x;
            
            err_back1 = norm(b - hb1, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            err_back2 = norm(b - hb2, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            ref_err_back = norm(b - rhb, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            err_bound = 2*(sqrt(2)+1)*sqrt(2^(depth + 1)+2^(depth-1))*eps;

            err_back_list1(i, j) = err_back_list1(i, j) + err_back1;
            err_back_list2(i, j) = err_back_list2(i, j) + err_back2;
            ref_err_back_list(i, j) = ref_err_back_list(i, j) + ref_err_back;
            err_bound_list(i, j) = err_bound_list(i, j) + err_bound;
        end
    end
end


save("results/prod3_err_back1.mat", 'err_back_list1');
save("results/prod3_err_back2.mat", 'err_back_list2');
save("results/prod3_ref_err_back.mat", 'ref_err_back_list');
save("results/prod3_bound.mat", 'err_bound_list');




%% kernel matrix 4
clear all
rng(0)

n_sample = 10;
x = unifrnd(-1, 1, 2000, 2);
y = unifrnd(-1, 1, 2000, 2);
kernel_mat = kernel4(x, y);

rng(42);
v =  unifrnd(-1, 1, n_sample, 2000);

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
ref_err_back_list  = zeros(n_eps, n_d);
err_bound_list = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, kernel_mat, depth, 10, 'svd', eps); 

        for k=1:n_sample
            x = v(k, :)';
            
            hb1 = mhdot(aphA, x, u2, 'dense');
            hb2 = mhdot(aphA, x, u4, 'dense');
            rhb = hdot(aphA, x, 'dense');
            
            b = kernel_mat * x;
            
            err_back1 = norm(b - hb1, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            err_back2 = norm(b - hb2, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            ref_err_back = norm(b - rhb, 'fro') / (norm(b, 'fro') * norm(kernel_mat, 'fro'));
            err_bound = 2*(sqrt(2)+1)*sqrt(2^(depth + 1)+2^(depth-1))*eps;

            err_back_list1(i, j) = err_back_list1(i, j) + err_back1;
            err_back_list2(i, j) = err_back_list2(i, j) + err_back2;
            ref_err_back_list(i, j) = ref_err_back_list(i, j) + ref_err_back;
            err_bound_list(i, j) = err_bound_list(i, j) + err_bound;
        end
    end
end


save("results/prod4_err_back1.mat", 'err_back_list1');
save("results/prod4_err_back2.mat", 'err_back_list2');
save("results/prod4_ref_err_back.mat", 'ref_err_back_list');
save("results/prod4_bound.mat", 'err_bound_list');

