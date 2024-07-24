%% kernel matrix 1
rng(0)

n_sample = 10;
x = rand(1, 2000);
y = rand(1, 2000);
kernel_d1 = kernel1(x, y);
v = rand(n_sample, 2000);

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

ref_err_forward_list  = zeros(n_eps, n_d);
ref_err_back_list  = zeros(n_eps, n_d);

err_forward_bound_list  = zeros(n_eps, n_d);
err_back_bound_list  = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, kernel_d1, depth, 10, 'svd', eps); 
        hA = hodlr(kernel_d1, depth, 10, 'svd', eps); 

        for k=1:n_sample
            x = v(k, :)';
            
            hb = hdot(aphA, x, 'dense');
            rhb = hdot(hA, x, 'dense');
            
            b = kernel_d1 * x;
            
            delta_H_bound = 2*(sqrt(2) + 1)*sqrt(2^(depth+1) + 2^(depth-1)) * eps * norm(kernel_d1, 'fro');
            
            err_forward = norm(b - hb, 'fro') / norm(b, 'fro');
            err_back = norm(b - hb, 'fro') / (norm(b, 'fro') * norm(kernel_d1, 'fro'));

            ref_err_forward = norm(b - rhb, 'fro') / norm(b, 'fro');
            ref_err_back = norm(b - rhb, 'fro') / (norm(b, 'fro') * norm(kernel_d1, 'fro'));
            
            err_forward_bound = delta_H_bound * norm(x, 'fro') / norm(b, 'fro');
            err_back_bound = delta_H_bound / norm(kernel_d1, 'fro');
    
            err_forward_list(i, j) = err_forward_list(i, j) + err_forward;
            err_back_list(i, j) = err_back_list(i, j) + err_back;
            
            ref_err_forward_list(i, j) = ref_err_forward_list(i, j) + ref_err_forward;
            ref_err_back_list(i, j) = ref_err_back_list(i, j) + ref_err_back;

            err_forward_bound_list(i, j) = err_forward_bound_list(i, j) + err_forward_bound;
            err_back_bound_list(i, j) = err_back_bound_list(i, j) + err_back_bound;
        end
    end
end


fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_forward_list(:, 1)./n_sample, ':c', 'Marker', '+', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,'--y', 'Marker', 'pentagram', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_list(:, 1)./n_sample, ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1)./n_sample,'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 1)./n_sample,':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 1)./n_sample,'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
yline(u1.u, '--')
yline(u2.u, '--')
yline(u3.u, '--')

h = legend('forward error (I)', ...
    'backward error (I)', ...
    'forward error (II)', ...
    'backward error (II)', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3)
legend boxoff
rect = [0.32, 0.835, .25, 0];
set(h, 'Position', rect)

t1 = text(8, u3.u, 'fp16')
t2 = text(8, u2.u, 'fp32')
t3 = text(8, u1.u, 'fp64')

t1.FontSize = fontSize;
t2.FontSize = fontSize;
t3.FontSize = fontSize;

[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel1_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_forward_list(:, 2)./n_sample, ':c', 'Marker', '+', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, ref_err_back_list(:, 2)./n_sample,'--y', 'Marker', 'pentagram', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_list(:, 2)./n_sample, ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 2)./n_sample,'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 2)./n_sample,':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 2)./n_sample,'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
yline(u1.u, '--')
yline(u2.u, '--')
yline(u3.u, '--')

h = legend('forward error (I)', ...
    'backward error (I)', ...
    'forward error (II)', ...
    'backward error (II)', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3)
legend boxoff
rect = [0.32, 0.835, .25, 0];
set(h, 'Position', rect)

t1 = text(8, u3.u, 'fp16')
t2 = text(8, u2.u, 'fp32')
t3 = text(8, u1.u, 'fp64')

t1.FontSize = fontSize;
t2.FontSize = fontSize;
t3.FontSize = fontSize;

[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel1_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_forward_list(:, 3)./n_sample, ':c', 'Marker', '+', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, ref_err_back_list(:, 3)./n_sample,'--y', 'Marker', 'pentagram', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_list(:, 3)./n_sample, ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 3)./n_sample,'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 3)./n_sample,':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 3)./n_sample,'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
yline(u1.u, '--')
yline(u2.u, '--')
yline(u3.u, '--')

h = legend('forward error (I)', ...
    'backward error (I)', ...
    'forward error (II)', ...
    'backward error (II)', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3)
legend boxoff
rect = [0.32, 0.835, .25, 0];
set(h, 'Position', rect)

t1 = text(8, u3.u, 'fp16')
t2 = text(8, u2.u, 'fp32')
t3 = text(8, u1.u, 'fp64')

t1.FontSize = fontSize;
t2.FontSize = fontSize;
t3.FontSize = fontSize;

[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel1_d8.pdf')
hold off




%% kernel matrix 2
clear all
rng(0)

n_sample = 10;
x = rand(1, 2000);
y = rand(1, 2000);
kernel_d1 = kernel2(x, y);
v = rand(n_sample, 2000);

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
        for k=1:n_sample
            x = v(k, :)';
            
            hb = hdot(aphA, x, 'dense');
            
            b = kernel_d1 * x;
            
            delta_H_bound = 2*(sqrt(2) + 1)*sqrt(2^(depth+1) + 2^(depth-1)) * eps * norm(kernel_d1, 'fro');
            
            err_forward = norm(b - hb, 'fro') / norm(b, 'fro');
            err_back = norm(b - hb, 'fro') / (norm(b, 'fro') * norm(kernel_d1, 'fro'));
            
            err_forward_bound = delta_H_bound * norm(x, 'fro') / norm(b, 'fro');
            err_back_bound = delta_H_bound / norm(kernel_d1, 'fro');
    
            err_forward_list(i, j) = err_forward_list(i, j) + err_forward;
            err_back_list(i, j) = err_back_list(i, j) + err_back;
    
            err_forward_bound_list(i, j) = err_forward_bound_list(i, j) + err_forward_bound;
            err_back_bound_list(i, j) = err_back_bound_list(i, j) + err_back_bound;
        end
    end
end


fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_forward_list(:, 1)./n_sample, ':c', 'Marker', '+', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,'--y', 'Marker', 'pentagram', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_list(:, 1)./n_sample, ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1)./n_sample,'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 1)./n_sample,':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 1)./n_sample,'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
yline(u1.u, '--')
yline(u2.u, '--')
yline(u3.u, '--')

h = legend('forward error (I)', ...
    'backward error (I)', ...
    'forward error (II)', ...
    'backward error (II)', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3)
legend boxoff
rect = [0.32, 0.835, .25, 0];
set(h, 'Position', rect)

t1 = text(8, u3.u, 'fp16')
t2 = text(8, u2.u, 'fp32')
t3 = text(8, u1.u, 'fp64')

t1.FontSize = fontSize;
t2.FontSize = fontSize;
t3.FontSize = fontSize;

[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel2_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_forward_list(:, 2)./n_sample, ':c', 'Marker', '+', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, ref_err_back_list(:, 2)./n_sample,'--y', 'Marker', 'pentagram', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_list(:, 2)./n_sample, ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 2)./n_sample,'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 2)./n_sample,':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 2)./n_sample,'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
yline(u1.u, '--')
yline(u2.u, '--')
yline(u3.u, '--')

h = legend('forward error (I)', ...
    'backward error (I)', ...
    'forward error (II)', ...
    'backward error (II)', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3)
legend boxoff
rect = [0.32, 0.835, .25, 0];
set(h, 'Position', rect)

t1 = text(8, u3.u, 'fp16')
t2 = text(8, u2.u, 'fp32')
t3 = text(8, u1.u, 'fp64')

t1.FontSize = fontSize;
t2.FontSize = fontSize;
t3.FontSize = fontSize;

[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel2_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_forward_list(:, 3)./n_sample, ':c', 'Marker', '+', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, ref_err_back_list(:, 3)./n_sample,'--y', 'Marker', 'pentagram', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_list(:, 3)./n_sample, ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 3)./n_sample,'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 3)./n_sample,':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 3)./n_sample,'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
yline(u1.u, '--')
yline(u2.u, '--')
yline(u3.u, '--')

h = legend('forward error (I)', ...
    'backward error (I)', ...
    'forward error (II)', ...
    'backward error (II)', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3)
legend boxoff
rect = [0.32, 0.835, .25, 0];
set(h, 'Position', rect)

t1 = text(8, u3.u, 'fp16')
t2 = text(8, u2.u, 'fp32')
t3 = text(8, u1.u, 'fp64')

t1.FontSize = fontSize;
t2.FontSize = fontSize;
t3.FontSize = fontSize;

[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel2_d8.pdf')
hold off


%% kernel matrix 3
clear all
rng(0)

n_sample = 10;
x = rand(2000, 3);
y = rand(2000, 3);
kernel_d1 = kernel3(x, y);
v = rand(n_sample, 2000);

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
        for k=1:n_sample
            x = v(k, :)';
            
            hb = hdot(aphA, x, 'dense');
            
            b = kernel_d1 * x;
            
            delta_H_bound = 2*(sqrt(2) + 1)*sqrt(2^(depth+1) + 2^(depth-1)) * eps * norm(kernel_d1, 'fro');
            
            err_forward = norm(b - hb, 'fro') / norm(b, 'fro');
            err_back = norm(b - hb, 'fro') / (norm(b, 'fro') * norm(kernel_d1, 'fro'));
            
            err_forward_bound = delta_H_bound * norm(x, 'fro') / norm(b, 'fro');
            err_back_bound = delta_H_bound / norm(kernel_d1, 'fro');
    
            err_forward_list(i, j) = err_forward_list(i, j) + err_forward;
            err_back_list(i, j) = err_back_list(i, j) + err_back;
    
            err_forward_bound_list(i, j) = err_forward_bound_list(i, j) + err_forward_bound;
            err_back_bound_list(i, j) = err_back_bound_list(i, j) + err_back_bound;
        end
    end
end


fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_forward_list(:, 1)./n_sample, ':c', 'Marker', '+', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,'--y', 'Marker', 'pentagram', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_list(:, 1)./n_sample, ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1)./n_sample,'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 1)./n_sample,':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 1)./n_sample,'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
yline(u1.u, '--')
yline(u2.u, '--')
yline(u3.u, '--')

h = legend('forward error (I)', ...
    'backward error (I)', ...
    'forward error (II)', ...
    'backward error (II)', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3)
legend boxoff
rect = [0.32, 0.835, .25, 0];
set(h, 'Position', rect)

t1 = text(8, u3.u, 'fp16')
t2 = text(8, u2.u, 'fp32')
t3 = text(8, u1.u, 'fp64')

t1.FontSize = fontSize;
t2.FontSize = fontSize;
t3.FontSize = fontSize;

[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel3_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_forward_list(:, 2)./n_sample, ':c', 'Marker', '+', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, ref_err_back_list(:, 2)./n_sample,'--y', 'Marker', 'pentagram', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_list(:, 2)./n_sample, ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 2)./n_sample,'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 2)./n_sample,':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 2)./n_sample,'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
yline(u1.u, '--')
yline(u2.u, '--')
yline(u3.u, '--')

h = legend('forward error (I)', ...
    'backward error (I)', ...
    'forward error (II)', ...
    'backward error (II)', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3)
legend boxoff
rect = [0.32, 0.835, .25, 0];
set(h, 'Position', rect)

t1 = text(8, u3.u, 'fp16')
t2 = text(8, u2.u, 'fp32')
t3 = text(8, u1.u, 'fp64')

t1.FontSize = fontSize;
t2.FontSize = fontSize;
t3.FontSize = fontSize;

[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel3_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_forward_list(:, 3)./n_sample, ':c', 'Marker', '+', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, ref_err_back_list(:, 3)./n_sample,'--y', 'Marker', 'pentagram', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_list(:, 3)./n_sample, ':r', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 3)./n_sample,'--g', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_forward_bound_list(:, 3)./n_sample,':k', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_bound_list(:, 3)./n_sample,'--b', 'Marker', 's',  'MarkerSize',10, 'Linewidth',2)
yline(u1.u, '--')
yline(u2.u, '--')
yline(u3.u, '--')

h = legend('forward error (I)', ...
    'backward error (I)', ...
    'forward error (II)', ...
    'backward error (II)', ...
    'forward error bound', 'backward error bound', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3)
legend boxoff
rect = [0.32, 0.835, .25, 0];
set(h, 'Position', rect)

t1 = text(8, u3.u, 'fp16')
t2 = text(8, u2.u, 'fp32')
t3 = text(8, u1.u, 'fp64')

t1.FontSize = fontSize;
t2.FontSize = fontSize;
t3.FontSize = fontSize;

[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel3_d8.pdf')
hold off