vareps = [1e-7, 1e-4, 1e-1]; % u2.u, u3.u
depths = [5, 8];


load("results/storage_eff1.mat");
load("results/storage_eff_ap1.mat");

fontSize = 19;

set(gcf, 'Position',  [10 10 400 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
ste = storage_eff./storage_eff_ap;
ba = bar(ste(2, :))
ba.FaceColor = 'flat';
ba.CData(1,:) = [.5 0 .5];
ba.CData(2,:) = [.5 1 .5];
ba.CData(3,:) = [.5 .5 1];

grid on;
ylim([0,2])
xticklabels(split(['$\varepsilon$=10^{-7}', ' $\varepsilon$=10^{-4}', ' $\varepsilon$=10^{-1}']))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
%legend(split(num2str(vareps,'%.e ')),'location','northeastoutside', BackgroundAlpha=.3)
exportgraphics(gca, 'figures/storage_P64.pdf')
hold off

%% ============================================

load("results/storage_eff2.mat");
load("results/storage_eff_ap2.mat");

fontSize = 19;

set(gcf, 'Position',  [10 10 400 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
ste = storage_eff./storage_eff_ap;
ba = bar(ste(2, :))
ba.FaceColor = 'flat';
ba.CData(1,:) = [.5 0 .5];
ba.CData(2,:) = [.5 1 .5];
ba.CData(3,:) = [.5 .5 1];

grid on;
ylim([0,2])
xticklabels(split(['$\varepsilon$=10^{-7}', ' $\varepsilon$=10^{-4}', ' $\varepsilon$=10^{-1}']))
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
%legend(split(num2str(vareps,'%.e ')),'location','northeastoutside', BackgroundAlpha=.3)
exportgraphics(gca, 'figures/storage_LeGresley_2508.pdf')
hold off

%% ============================================

load("results/storage_eff3.mat");
load("results/storage_eff_ap3.mat");

fontSize = 19;

set(gcf, 'Position',  [10 10 400 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
ste = storage_eff./storage_eff_ap;
ba = bar(ste(2, :))
ba.FaceColor = 'flat';
ba.CData(1,:) = [.5 0 .5];
ba.CData(2,:) = [.5 1 .5];
ba.CData(3,:) = [.5 .5 1];

grid on;
ylim([0,2])
xticklabels(split(['$\varepsilon$=10^{-7}', ' $\varepsilon$=10^{-4}', ' $\varepsilon$=10^{-1}']))

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
%legend(split(num2str(vareps,'%.e ')),'location','northeastoutside', BackgroundAlpha=.3)
exportgraphics(gca, 'figures/storage_ex37.pdf')
hold off

%% ============================================

load("results/storage_eff4.mat");
load("results/storage_eff_ap4.mat");

fontSize = 19;

set(gcf, 'Position',  [10 10 400 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
ste = storage_eff./storage_eff_ap;
ba = bar(ste(2, :))
ba.FaceColor = 'flat';
ba.CData(1,:) = [.5 0 .5];
ba.CData(2,:) = [.5 1 .5];
ba.CData(3,:) = [.5 .5 1];

grid on;
ylim([0,2])
xticklabels(split(['$\varepsilon$=10^{-7}', ' $\varepsilon$=10^{-4}', ' $\varepsilon$=10^{-1}']))

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
%legend(split(num2str(vareps,'%.e ')),'location','northeastoutside', BackgroundAlpha=.3)
exportgraphics(gca, 'figures/storage_1138_bus.pdf')
hold off

%% ============================================

load("results/storage_eff5.mat");
load("results/storage_eff_ap5.mat");

fontSize = 19;

set(gcf, 'Position',  [10 10 400 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
ste = storage_eff./storage_eff_ap;
ba = bar(ste(2, :))
ba.FaceColor = 'flat';
ba.CData(1,:) = [.5 0 .5];
ba.CData(2,:) = [.5 1 .5];
ba.CData(3,:) = [.5 .5 1];

grid on;
ylim([0,2])
xticklabels(split(['$\varepsilon$=10^{-7}', ' $\varepsilon$=10^{-4}', ' $\varepsilon$=10^{-1}']))

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
%legend(split(num2str(vareps,'%.e ')),'location','northeastoutside', BackgroundAlpha=.3)
exportgraphics(gca, 'figures/storage_cavity18.pdf')
hold off

%% ============================================
load("results/storage_eff6.mat");
load("results/storage_eff_ap6.mat");

fontSize = 19;

set(gcf, 'Position',  [10 10 400 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
ste = storage_eff./storage_eff_ap;
ba = bar(ste(2, :))
ba.FaceColor = 'flat';
ba.CData(1,:) = [.5 0 .5];
ba.CData(2,:) = [.5 1 .5];
ba.CData(3,:) = [.5 .5 1];

grid on;
ylim([0,2])
xticklabels(split(['$\varepsilon$=10^{-7}', ' $\varepsilon$=10^{-4}', ' $\varepsilon$=10^{-1}']))

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
%legend(split(num2str(vareps,'%.e ')),'location','northeastoutside', BackgroundAlpha=.3)
exportgraphics(gca, 'figures/storage_psmigr_1.pdf')
hold off


%% ============================================
load("results/storage_eff7.mat");
load("results/storage_eff_ap7.mat");

fontSize = 19;

set(gcf, 'Position',  [10 10 400 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
ste = storage_eff./storage_eff_ap;
ba = bar(ste(2, :))
ba.FaceColor = 'flat';
ba.CData(1,:) = [.5 0 .5];
ba.CData(2,:) = [.5 1 .5];
ba.CData(3,:) = [.5 .5 1];

grid on;
ylim([0,2])
xticklabels(split(['$\varepsilon$=10^{-7}', ' $\varepsilon$=10^{-4}', ' $\varepsilon$=10^{-1}']))

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
%legend(split(num2str(vareps,'%.e ')),'location','northeastoutside', BackgroundAlpha=.3)
exportgraphics(gca, 'figures/storage_breasttissue_10NN.pdf')
hold off


%% ============================================
load("results/storage_eff8.mat");
load("results/storage_eff_ap8.mat");

fontSize = 19;

set(gcf, 'Position',  [10 10 400 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
ste = storage_eff./storage_eff_ap;
ba = bar(ste(2, :))
ba.FaceColor = 'flat';
ba.CData(1,:) = [.5 0 .5];
ba.CData(2,:) = [.5 1 .5];
ba.CData(3,:) = [.5 .5 1];

grid on;
ylim([0,2])
xticklabels(split(['$\varepsilon$=10^{-7}', ' $\varepsilon$=10^{-4}', ' $\varepsilon$=10^{-1}']))

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
%legend(split(num2str(vareps,'%.e ')),'location','northeastoutside', BackgroundAlpha=.3)
exportgraphics(gca, 'figures/storage_saylr3.pdf')
hold off



%% ============================================
load("results/storage_eff9.mat");
load("results/storage_eff_ap9.mat");

fontSize = 19;

set(gcf, 'Position',  [10 10 400 500])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
ste = storage_eff./storage_eff_ap;
ba = bar(ste(2, :))
ba.FaceColor = 'flat';
ba.CData(1,:) = [.5 0 .5];
ba.CData(2,:) = [.5 1 .5];
ba.CData(3,:) = [.5 .5 1];

grid on;
ylim([0,2])
xticklabels(split(['$\varepsilon$=10^{-7}', ' $\varepsilon$=10^{-4}', ' $\varepsilon$=10^{-1}']))

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',fontSize) % ,'FontWeight','bold'
%legend(split(num2str(vareps,'%.e ')),'location','northeastoutside', BackgroundAlpha=.3)
exportgraphics(gca, 'figures/storage_bcsstk08.pdf')
hold off