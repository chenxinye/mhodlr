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

err_back_list  = zeros(n_eps, n_d);
ref_err_back_list  = zeros(n_eps, n_d);


for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        [L1, U1] = hlu(aphA, 'hodlr');
        [L2, U2] = mhlu(aphA, u2, 'hodlr');
        
        b = A * x;
        
        err_back1 = norm(b - hb1, 'fro') / (norm(b, 'fro') * norm(A, 'fro'));
        ref_err_back = norm(b - rhb, 'fro') / (norm(b, 'fro') * norm(A, 'fro'));

        err_back_list(i, j) = err_back_list(i, j) + err_back1;
        ref_err_back_list(i, j) = ref_err_back_list(i, j) + ref_err_back;

    end
end


fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

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

semilogy(1:n_eps, ref_err_back_list(:, 2)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 2)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel1_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel1_d8.pdf')
hold off




%% kernel matrix 2
rng(0)

A = load('data/LeGresley_2508.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

depths = [2, 5, 8];
vareps = [1e-14, 1e-12, 1e-10, 1e-08, 1e-06, 1e-04, 1e-02, 1e-00];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_back_list  = zeros(n_eps, n_d);
ref_err_back_list  = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        for k=1:n_sample
            x = v(k, :)';
            
            hb1 = mhdot(aphA, x, u2, 'dense');
            rhb = hdot(aphA, x, 'dense');
            
            b = A * x;
            
            err_back1 = norm(b - hb1, 'fro') / (norm(b, 'fro') * norm(A, 'fro'));
            ref_err_back = norm(b - rhb, 'fro') / (norm(b, 'fro') * norm(A, 'fro'));
    
            err_back_list(i, j) = err_back_list(i, j) + err_back1;
            ref_err_back_list(i, j) = ref_err_back_list(i, j) + ref_err_back;
        end
    end
end


fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')));

lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel2_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 2)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 2)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel2_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel2_d8.pdf')
hold off


%% kernel matrix 3
clear all
rng(0)

A = load('data/breasttissue_10NN.mat');
A =  A.Problem.A;

u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

depths = [2, 5, 8];
vareps = [1e-14, 1e-12, 1e-10, 1e-08, 1e-06, 1e-04, 1e-02, 1e-00];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

err_back_list  = zeros(n_eps, n_d);
ref_err_back_list  = zeros(n_eps, n_d);

for i=1:n_eps
    for j=1:n_d
        eps = vareps(i);
        depth = depths(j);
        
        aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
        
        hb1 = mhdot(aphA, x, u2, 'dense');
        rhb = hdot(aphA, x, 'dense');
        
        b = A * x;
        
        err_back1 = norm(b - hb1, 'fro') / (norm(b, 'fro') * norm(A, 'fro'));
        ref_err_back = norm(b - rhb, 'fro') / (norm(b, 'fro') * norm(A, 'fro'));

        err_back_list(i, j) = err_back_list(i, j) + err_back1;
        ref_err_back_list(i, j) = ref_err_back_list(i, j) + ref_err_back;
    
    end
end


fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel3_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 2)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 2)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel3_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,':g', 'Marker', 'o', 'MarkerSize',10, 'Linewidth',2)
hold on 
semilogy(1:n_eps, err_back_list(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',10, 'Linewidth',2)
hold on 

h = legend('backward error (fp64)', ...
           'backward error (fp32)', ...
     'NumColumns',2, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.36, 0.835, .25, 0];
set(h, 'Position', rect);

[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')));


lx = xlabel('$\varepsilon$');
set(lx,'interpreter','latex');
lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/matvecprod_kernel3_d8.pdf')
hold off