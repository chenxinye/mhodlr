%% kernel matrix 1
rng(0)

x = rand(1, 2000);
y = rand(1, 2000);
kernel_d1 = kernel1(x, y);
v = rand(1, 2000);

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

depths = [2, 5, 8];
vareps = [10e-9, 10e-08, 10e-07, 10e-06, 10e-5, 10e-04, 10e-03, 10e-02];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_forward_list  = zeros(n_eps, n_d);
err_back_list  = zeros(n_eps, n_d);

err_forward_bound_list  = zeros(n_eps, n_d);
err_back_bound_list  = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);

        aphA = amphodlr(u_chain, kernel_d1, depth, 10, 'svd', eps); 
        hb = hdot(aphA, v', 'dense');
        
        b = kernel_d1 * v';
        
        delta_H_bound = 2*(sqrt(2) + 1)*sqrt(2^(depth+1) + 2^(depth-1)) * eps * norm(kernel_d1, 'fro');
        
        err_forward = norm(b - hb, 'fro') / norm(b, 'fro');
        err_back = norm(b - hb, 'fro') / (norm(b, 'fro') * norm(kernel_d1, 'fro'));
        
        err_forward_bound = delta_H_bound * norm(v, 'fro') / norm(b, 'fro');
        err_back_bound = delta_H_bound / norm(kernel_d1, 'fro');

        err_forward_list(i, j) = err_forward;
        err_back_list(i, j) = err_back;

        err_forward_bound_list(i, j) = err_forward_bound;
        err_back_bound_list(i, j) = err_back_bound;
    end
end


fontSize = 18;
set(gcf, 'Position',  [10 10 800 650])
semilogy(1:n_eps, err_forward_list(:, 1), ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1),'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 1),':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 1),'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)

legend('forward error', ...
    'backward error', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize)

[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel1_d2.pdf')
hold off


set(gcf, 'Position',  [10 10 800 650])
semilogy(1:n_eps, err_forward_list(:, 2), ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 2),'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 2),':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 2),'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)

legend('forward error', ...
    'backward error', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize)

[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel1_d5.pdf')
hold off


set(gcf, 'Position',  [10 10 800 650])
semilogy(1:n_eps, err_forward_list(:, 3), ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 3),'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 3),':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 3),'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)

legend('forward error', ...
    'backward error', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize)

[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel1_d8.pdf')
hold off




%% kernel matrix 2
rng(0)

x = rand(1, 2000);
y = rand(1, 2000);
kernel_d1 = kernel2(x, y);
v = rand(1, 2000);

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

depths = [2, 5, 8];
vareps = [10e-9, 10e-08, 10e-07, 10e-06, 10e-5, 10e-04, 10e-03, 10e-02];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_forward_list  = zeros(n_eps, n_d);
err_back_list  = zeros(n_eps, n_d);

err_forward_bound_list  = zeros(n_eps, n_d);
err_back_bound_list  = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);

        aphA = amphodlr(u_chain, kernel_d1, depth, 10, 'svd', eps); 
        hb = hdot(aphA, v', 'dense');
        
        b = kernel_d1 * v';
        
        delta_H_bound = 2*(sqrt(2) + 1)*sqrt(2^(depth+1) + 2^(depth-1)) * eps * norm(kernel_d1, 'fro');
        
        err_forward = norm(b - hb, 'fro') / norm(b, 'fro');
        err_back = norm(b - hb, 'fro') / (norm(b, 'fro') * norm(kernel_d1, 'fro'));
        
        err_forward_bound = delta_H_bound * norm(v, 'fro') / norm(b, 'fro');
        err_back_bound = delta_H_bound / norm(kernel_d1, 'fro');

        err_forward_list(i, j) = err_forward;
        err_back_list(i, j) = err_back;

        err_forward_bound_list(i, j) = err_forward_bound;
        err_back_bound_list(i, j) = err_back_bound;
    end
end

fontSize = 18;
set(gcf, 'Position',  [10 10 800 650])
semilogy(1:n_eps, err_forward_list(:, 1), ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1),'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 1),':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 1),'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)

legend('forward error', ...
    'backward error', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize)

[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel2_d2.pdf')
hold off


set(gcf, 'Position',  [10 10 800 650])
semilogy(1:n_eps, err_forward_list(:, 2), ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 2),'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 2),':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 2),'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)


legend('forward error', ...
    'backward error', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize)

[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel2_d5.pdf')
hold off


set(gcf, 'Position',  [10 10 800 650])
semilogy(1:n_eps, err_forward_list(:, 3), ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 3),'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 3),':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 3),'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)


legend('forward error', ...
    'backward error', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize)

[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel2_d8.pdf')
hold off