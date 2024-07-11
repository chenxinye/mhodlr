%% P64_cs128
vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/blrLU/root_P64_cs128.mat');
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
    for j=1:n
        d = depths(i);
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        
        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

fontSize = 13;

set(gcf, 'Position',  [10 10 600 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
bar(storage_eff./storage_eff_ap)
xticklabels(split(['$\ell$=5', ' $\ell$=8']))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
legend(split(num2str(vareps,'%.e ')),'location','northeastoutside')
exportgraphics(gca, 'figures/P64_cs128_storage.pdf')
hold off


%% LeGresley_2508

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
    for j=1:n
        d = depths(i);
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        
        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

fontSize = 13;

set(gcf, 'Position',  [10 10 600 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
bar(storage_eff./storage_eff_ap)
xticklabels(split(['$\ell$=5', ' $\ell$=8']))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
legend(split(num2str(vareps,'%.e ')),'location','northeastoutside')
exportgraphics(gca, 'figures/LeGresley_2508_storage.pdf')
hold off


%% breasttissue_10NN

vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/breasttissue_10NN.mat');
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
    for j=1:n
        d = depths(i);
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        
        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

fontSize = 13;

set(gcf, 'Position',  [10 10 600 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
bar(storage_eff./storage_eff_ap)
xticklabels(split(['$\ell$=5', ' $\ell$=8']))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
legend(split(num2str(vareps,'%.e ')),'location','northeastoutside')
exportgraphics(gca, 'figures/breasttissue_10NN_storage.pdf')
hold off



%% 1138_bus

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
    for j=1:n
        d = depths(i);
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        
        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

fontSize = 13;

set(gcf, 'Position',  [10 10 600 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
bar(storage_eff./storage_eff_ap)
xticklabels(split(['$\ell$=5', ' $\ell$=8']))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
legend(split(num2str(vareps,'%.e ')),'location','northeastoutside')
exportgraphics(gca, 'figures/1138_bus_storage.pdf')
hold off



%% lshp1009

vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/lshp1009.mat');
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
    for j=1:n
        d = depths(i);
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        
        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

fontSize = 13;

set(gcf, 'Position',  [10 10 600 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
bar(storage_eff./storage_eff_ap)
xticklabels(split(['$\ell$=5', ' $\ell$=8']))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
legend(split(num2str(vareps,'%.e ')),'location','northeastoutside')
exportgraphics(gca, 'figures/lshp1009_storage.pdf')
hold off


%% west2021.mat

vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];

A = load('data/west2021.mat');
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
    for j=1:n
        d = depths(i);
        eps = vareps(j);
        hA = hodlr(A, d, 2, 'svd', eps);
        rA = recover(hA);
        
        aphA = amphodlr(u_chain, A, d, 2, 'svd', eps);
        aprA = recover(aphA);
        
        [n1, n2] = size(A);
        storage_eff(i, j) = hstorage(hA) / (n1*n2*64);
        storage_eff_ap(i, j) = hstorage(aphA) / (n1*n2*64);
        disp([storage_eff(i, j), storage_eff_ap(i, j)])
    end
end

fontSize = 13;

set(gcf, 'Position',  [10 10 600 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
bar(storage_eff./storage_eff_ap)
xticklabels(split(['$\ell$=5', ' $\ell$=8']))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
legend(split(num2str(vareps,'%.e ')),'location','northeastoutside')
exportgraphics(gca, 'figures/west2021_storage.pdf')
hold off