clear all 
delete(gca)
delete(gcf)

vareps = [10e-9, 10e-08, 10e-07, 10e-06, 10e-5, 10e-04, 10e-03, 10e-02]; % u2.u, u3.u
depths = [2,  5,  8];

n_d = size(depths, 2);
n_eps = size(vareps, 2);

%% P64
load("results/rcerr_err_list1.mat");
load("results/rcerr_err_list_amp1.mat");
load("results/rcerr_err_bound_list1.mat");

rect = [0.5, 0, 0, 0.05];
fontSize = 18;
lineWidth = 2;
markerSize = 15;

set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 1), '-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_list_amp(:, 1),'--r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 1),':k', 'Marker', '.',  'MarkerSize',markerSize, 'Linewidth',lineWidth)

legend('fp64', ...
    'adaptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
legend boxoff
[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
ylim([1e-10, 1e-0]);
yticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-0]);
grid on;
exportgraphics(gca, 'figures/P64_depth1.pdf')
hold off


set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 2), '-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_list_amp(:, 2),'--r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 2),':k', 'Marker', '.',  'MarkerSize',markerSize, 'Linewidth',lineWidth)

legend('fp64', ...
    'adaptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
legend boxoff
[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
ylim([1e-10, 1e-0]);
yticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-0]);
grid on;
exportgraphics(gca, 'figures/P64_depth2.pdf')
hold off

set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 3), '-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_list_amp(:, 3),'--r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 3),':k', 'Marker', '.',  'MarkerSize',markerSize, 'Linewidth',lineWidth)

legend('fp64', ...
    'adaptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
legend boxoff
[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
ylim([1e-10, 1e-0]);
yticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-0]);
grid on;
exportgraphics(gca, 'figures/P64_depth3.pdf')
hold off
delete(gcf)
delete(gca)


% ex37
load("results/rcerr_err_list2.mat");
load("results/rcerr_err_list_amp2.mat");
load("results/rcerr_err_bound_list2.mat");


rect = [0.5, 0, 0, 0.05];
fontSize = 18;
lineWidth = 2;
markerSize = 15;

set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 1), '-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_list_amp(:, 1),'--r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 1),':k', 'Marker', '.',  'MarkerSize',markerSize, 'Linewidth',lineWidth)

legend('fp64', ...
    'adaptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
legend boxoff
[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/ex37_depth1.pdf')
hold off


set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 2), '-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_list_amp(:, 2),'--r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 2),':k', 'Marker', '.',  'MarkerSize',markerSize, 'Linewidth',lineWidth)

legend('fp64', ...
    'adaptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
legend boxoff
[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/ex37_depth2.pdf')
hold off


set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 3), '-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_list_amp(:, 3),'--r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 3),':k', 'Marker', '.',  'MarkerSize',markerSize, 'Linewidth',lineWidth)

legend('fp64', ...
    'adaptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
legend boxoff
[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/ex37_depth3.pdf')
hold off
delete(gcf)
delete(gca)





% 1138_bus
load("results/rcerr_err_list3.mat");
load("results/rcerr_err_list_amp3.mat");
load("results/rcerr_err_bound_list3.mat");

rect = [0.5, 0, 0, 0.05];
fontSize = 18;
lineWidth = 2;
markerSize = 15;

set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 1), '-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_list_amp(:, 1),'--r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 1),':k', 'Marker', '.',  'MarkerSize',markerSize, 'Linewidth',lineWidth)

legend('fp64', ...
    'adaptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
legend boxoff
[l, s] = title('$\ell$=2');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/1138_bus_depth1.pdf')
hold off

set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 2), '-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_list_amp(:, 2),'--r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 2),':k', 'Marker', '.',  'MarkerSize',markerSize, 'Linewidth',lineWidth)

legend('fp64', ...
    'adaptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
legend boxoff
[l, s] = title('$\ell$=5');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/1138_bus_depth2.pdf')
hold off

set(gcf, 'Position',  [10 10 800 600])
semilogy(1:n_eps, err_list(:, 3), '-b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_list_amp(:, 3),'--r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_bound_list(:, 3),':k', 'Marker', '.',  'MarkerSize',markerSize, 'Linewidth',lineWidth)

legend('fp64', ...
    'adaptive precision', ...
    'error bound', 'Position', rect, 'NumColumns',3, 'FontSize', fontSize)
legend boxoff
[l, s] = title('$\ell$=8');
xticklabels(split(num2str(vareps,'%.e ')))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
set(l,'interpreter','latex');
l.FontSize = fontSize+12;
grid on;
exportgraphics(gca, 'figures/1138_bus_depth3.pdf')
hold off