%% kernel matrix 1
rng(0)

A = load('data/blrLU/root_P64_cs128.mat');
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
        [L3, U3] = mhlu(aphA, u3, 'hodlr');

        recover_LU1 = hdot(L1, U1, 'dense');
        recover_LU2 = hdot(L2, U2, 'dense');
        recover_LU3 = hdot(L3, U3, 'dense');

        err_back1 = norm(A - recover_LU1, 'fro') / norm(A, 'fro'));
        err_back2 = norm(A - recover_LU2, 'fro') / norm(A, 'fro'));
        err_back3 = norm(A - recover_LU3, 'fro') / norm(A, 'fro'));

        err_back_list1(i, j) = err_back1;
        err_back_list2(i, j) = err_back2;
        err_back_list3(i, j) = err_back3;
    end
end


fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, err_back_list1(:, 1)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list2(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu1_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, err_back_list1(:, 2)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list2(:, 2)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu1_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, err_back_list1(:, 1)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list2(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu1_d8.pdf')
hold off

