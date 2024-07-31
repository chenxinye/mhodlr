
n_sample = 10;
depths = [2, 5, 8];
vareps = [1e-14, 1e-12, 1e-10, 1e-08, 1e-06, 1e-04, 1e-02];
lineWidth = 2;
markerSize = 15;

n_d = size(depths, 2);
n_eps = size(vareps, 2);

%% kernel matrix 1
fontSize = 19;
load("results/prod1_err_back1.mat");
load("results/prod1_err_back2.mat");
load("results/prod1_ref_err_back.mat");
load("results/prod2_bound.mat");

%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 1)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 1)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 


h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-22, 1e-0+0.5]);
yticks([1e-22, 1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');_
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel1_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 2)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 2)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 2)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 2)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-22, 1e-0+0.5]);
yticks([1e-22, 1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel1_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 3)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 3)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 3)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 3)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-22, 1e-0+0.5]);
yticks([1e-22, 1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel1_d8.pdf')
hold off



%% kernel matrix 2
fontSize = 19;
load("results/prod2_err_back1.mat");
load("results/prod2_err_back2.mat");
load("results/prod2_ref_err_back.mat");
load("results/prod2_bound.mat");

%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 1)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 1)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-18, 1e-0+0.5]);
yticks([1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel2_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 2)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 2)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 2)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 2)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-18, 1e-0+0.5]);
yticks([1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel2_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 3)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 3)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 3)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 3)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-18, 1e-0+0.5]);
yticks([1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel2_d8.pdf')
hold off




%% kernel matrix 3
fontSize = 19;
load("results/prod3_err_back1.mat");
load("results/prod3_err_back2.mat");
load("results/prod3_ref_err_back.mat");
load("results/prod3_bound.mat");


%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 1)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 1)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-18, 1e-0+0.5]);
yticks([1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel3_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 2)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 2)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 2)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 2)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-18, 1e-0+0.5]);
yticks([1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel3_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 3)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 3)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 3)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 3)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-18, 1e-0+0.5]);
yticks([1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel3_d8.pdf')
hold off



%% kernel matrix 4
fontSize = 19;
load("results/prod4_err_back1.mat");
load("results/prod4_err_back2.mat");
load("results/prod4_ref_err_back.mat");
load("results/prod4_bound.mat");



%% plot 1
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 1)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 1)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 1)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 1)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-18, 1e-0+0.5]);
yticks([1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel4_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 2)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 2)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 2)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 2)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-18, 1e-0+0.5]);
yticks([1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel4_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650])

semilogy(1:n_eps, ref_err_back_list(:, 3)./n_sample,'-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list1(:, 3)./n_sample,'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 3)./n_sample,'--g', 'Marker', 's', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 3)./n_sample,':k', 'Marker', '.', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bf16', ...
           'error bound', ...
     'NumColumns',4, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.5, 0, 0, 0.05];
set(h, 'Position', rect);

%[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-18, 1e-0+0.5]);
yticks([1e-18, 1e-12, 1e-6, 1e-0]);

%lx = xlabel('$\varepsilon$');
%set(lx,'interpreter','latex');
%lx.FontSize = fontSize;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
% set(l,'interpreter','latex');
% l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/matvecprod_kernel4_d8.pdf')
hold off