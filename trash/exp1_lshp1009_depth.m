%% global error

A = load('data/lshp1009.mat');
%LeGresley_2508.mat
%LeGresley_4908.mat
A =  A.Problem.A;
disp(size(A));
%% 'data/LeGresley_2508.mat'
%% mix precision
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

rect = [0.5, 0, 0, 0.05];
fontSize = 15;
set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 1), '-b', 'Marker', 'o', 'MarkerSize',8, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_list_amp(:, 1),'--r', 'Marker', '*', 'MarkerSize', 8, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_bound_list(:, 1),':k', 'Marker', 's',  'MarkerSize',8, 'Linewidth',2)

legend('fp64', ...
    'adptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+7;
exportgraphics(gca, 'figures/lshp1009_depth=1.pdf')
hold off

set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 2), '-b', 'Marker', 'o', 'MarkerSize',8, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_list_amp(:, 2),'--r', 'Marker', '*', 'MarkerSize', 8, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_bound_list(:, 2),':k', 'Marker', 's',  'MarkerSize',8, 'Linewidth',2)

legend('fp64', ...
    'adptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+7;
exportgraphics(gca, 'figures/lshp1009_depth=2.pdf')
hold off


set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 3), '-b', 'Marker', 'o', 'MarkerSize',8, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_list_amp(:, 3),'--r', 'Marker', '*', 'MarkerSize', 8, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_bound_list(:, 3),':k', 'Marker', 's',  'MarkerSize',8, 'Linewidth',2)

legend('fp64', ...
    'adptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+7;
exportgraphics(gca, 'figures/lshp1009_depth=3.pdf')
hold off