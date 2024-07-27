depths = [2, 5, 8];
vareps = [1e-14, 1e-12, 1e-10, 1e-08, 1e-06, 1e-04, 1e-02];

n_d = size(depths, 2);
n_eps = size(vareps, 2);
lineWidth = 1.5;
markerSize = 13;

%% matrix 1
load("results/lu1_depth2.mat");
load("results/lu1_depth5.mat");
load("results/lu1_depth8.mat");

fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 1),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 1),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 1),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns',3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=2');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-15,  1])
yticks([1e-15, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu1_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 2),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 2),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 2),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
     'NumColumns', 3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=5');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-15,  1])
yticks([1e-15, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu1_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 3),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 3),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 3),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns',3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=8');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-15,  1])
yticks([1e-15, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu1_d8.pdf')
hold off


% ---------------------------------------------------------------------------------
%% matrix 2
load("results/lu2_depth2.mat");
load("results/lu2_depth5.mat");
load("results/lu2_depth8.mat");

fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 1),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 1),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 1),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns',3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=2');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-15,  1])
yticks([1e-15, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu2_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 2),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 2),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 2),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns', 3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=5');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-15,  1])
yticks([1e-15, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu2_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 3),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 3),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 3),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns',3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=8');
set(l,'interpreter','latex');
xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-15,  1])
yticks([1e-15, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu2_d8.pdf')
hold off


% -----------------------------------------------------
%% matrix 3
load("results/lu3_depth2.mat");
load("results/lu3_depth5.mat");
load("results/lu3_depth8.mat");


fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 1),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 1),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 1),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns',3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=2');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-16,  1])
yticks([1e-16, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu3_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 2),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 2),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 2),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns', 3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=5');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-16,  1])
yticks([1e-16, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu3_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 3),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 3),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 3),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns',3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=8');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-16,  1])
yticks([1e-16, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu3_d8.pdf')
hold off



% ---------------------------------------------------------------------------------
%% matrix 4
load("results/lu4_depth2.mat");
load("results/lu4_depth5.mat");
load("results/lu4_depth8.mat");


fontSize = 18;
%% plot 1
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 1),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 1),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 1),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns',3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=2');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-15,  1])
yticks([1e-15, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu4_d2.pdf')
hold off

%% plot 2
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 2),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 2),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 2),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns', 3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=5');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-15,  1])
yticks([1e-15, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu4_d5.pdf')
hold off

%% plot 3
set(gcf, 'Position',  [10 10 800 650]);

semilogy(1:n_eps, err_back_list1(:, 3),':b', 'Marker', 'o', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list2(:, 3),'-.r', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 
semilogy(1:n_eps, err_back_list3(:, 3),'-.g', 'Marker', '*', 'MarkerSize',markerSize, 'Linewidth',lineWidth)
hold on 

h = legend('fp64', ...
           'fp32', ...
           'bfp16', ...
     'NumColumns',3, 'Location', 'Best', 'FontSize', fontSize, BackgroundAlpha=.3);
legend boxoff
rect = [0.38, 0.875, .25, 0];
set(h, 'Position', rect);

[l, ~] = title('$\ell$=8');
set(l,'interpreter','latex');

xticklabels(split(num2str(vareps,'%.e ')));
ylim([1e-15,  1])
yticks([1e-15, 1e-10, 1e-5, 1e-0]);
%extraInputs = {'interpreter','latex','fontsize',fontSize};
%lx = xlabel('$\varepsilon$', extraInputs{:});

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize); % ,'FontWeight','bold'
l.FontSize = fontSize+12;
exportgraphics(gca, 'figures/lu4_d8.pdf')
hold off
