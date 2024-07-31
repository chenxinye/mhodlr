clear all
%% P64_cs128
disp('------------P64_cs128------------')
vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/3-5000/root_P64_cs128.mat');
A =  A.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

m = size(depths, 2);
n = size(vareps, 2);
storage_eff = zeros(m, n);
storage_eff_ap = zeros(m, n);


for i=1:m
    d = depths(i);
    for j=1:n
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        disp(aphA.precIndex)

        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

save("results/storage_eff1.mat", 'storage_eff');
save("results/storage_eff_ap1.mat", 'storage_eff_ap');


%% LeGresley_2508
disp('------------LeGresley_2508------------')
vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/LeGresley_2508.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

m = size(depths, 2);
n = size(vareps, 2);
storage_eff = zeros(m, n);
storage_eff_ap = zeros(m, n);

for i=1:m
    d = depths(i);
    for j=1:n
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        disp(aphA.precIndex)
        
        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

save("results/storage_eff2.mat", 'storage_eff');
save("results/storage_eff_ap2.mat", 'storage_eff_ap');


%% ex37
disp('------------ex37------------')
vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/3-5000/ex37.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

m = size(depths, 2);
n = size(vareps, 2);
storage_eff = zeros(m, n);
storage_eff_ap = zeros(m, n);

for i=1:m
    d = depths(i);
    for j=1:n
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        disp(aphA.precIndex)

        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

save("results/storage_eff3.mat", 'storage_eff');
save("results/storage_eff_ap3.mat", 'storage_eff_ap');





%% 1138_bus
disp('------------1138_bus------------')
vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/1138_bus.mat');
A =  A.Problem.A;


u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

m = size(depths, 2);
n = size(vareps, 2);
storage_eff = zeros(m, n);
storage_eff_ap = zeros(m, n);

for i=1:m
    d = depths(i);
    for j=1:n
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        disp(aphA.precIndex)

        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

save("results/storage_eff4.mat", 'storage_eff');
save("results/storage_eff_ap4.mat", 'storage_eff_ap');




%% cavity18
disp('------------cavity18------------')
vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/3-5000/cavity18.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

m = size(depths, 2);
n = size(vareps, 2);
storage_eff = zeros(m, n);
storage_eff_ap = zeros(m, n);

for i=1:m
    d = depths(i);
    for j=1:n
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        disp(aphA.precIndex)

        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

save("results/storage_eff5.mat", 'storage_eff');
save("results/storage_eff_ap5.mat", 'storage_eff_ap');


%% psmigr_1
disp('------------psmigr_1------------')
vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/3-5000/psmigr_1.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

m = size(depths, 2);
n = size(vareps, 2);
storage_eff = zeros(m, n);
storage_eff_ap = zeros(m, n);

for i=1:m
    d = depths(i);
    for j=1:n
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        disp(aphA.precIndex)

        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

save("results/storage_eff6.mat", 'storage_eff');
save("results/storage_eff_ap6.mat", 'storage_eff_ap');



%% saylr3
disp('------------saylr3------------')
vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/saylr3.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

m = size(depths, 2);
n = size(vareps, 2);
storage_eff = zeros(m, n);
storage_eff_ap = zeros(m, n);

for i=1:m
    d = depths(i);
    for j=1:n
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        disp(aphA.precIndex)

        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

save("results/storage_eff8.mat", 'storage_eff');
save("results/storage_eff_ap8.mat", 'storage_eff_ap');




%% bcsstk08
disp('------------bcsstk08------------')
vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/bcsstk08.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

m = size(depths, 2);
n = size(vareps, 2);
storage_eff = zeros(m, n);
storage_eff_ap = zeros(m, n);

for i=1:m
    d = depths(i);
    for j=1:n
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        disp(aphA.precIndex)

        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

save("results/storage_eff9.mat", 'storage_eff');
save("results/storage_eff_ap9.mat", 'storage_eff_ap');

