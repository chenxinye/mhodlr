clear all
A = load('data/saylr3.mat');
%LeGresley_2508.mat
%LeGresley_4908.mat
A =  A.Problem.A;
%% 'data/LeGresley_2508.mat'
%% mix precision
u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

u_chain = prec_chain(u1, u2, u3, u4, u5);

epsilon = 1e-4; % 1e-1
depth = 6;
aphA = amphodlr(u_chain, A, depth, 10, 'svd', epsilon); 
aprA = recover(aphA);

disp(aphA);
disp(epsilon);

norm(aprA - A, 'fro') / norm(A, 'fro')
bound_err = (2 * sqrt(2 * aphA.bottom_level) + 1) * epsilon 

pA = compute_hmat_prec(aphA);
h = heatmap(pA,'CellLabelColor','none');

h.ColorbarVisible = 'off';
h.GridVisible = 'off';
h.Colormap = spring;
set(gcf, 'Position',  [10 10 700 600])
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
exportgraphics(gca,'figures/precssaylr3.pdf')

nA = compute_hmat_norm(aphA);
h = heatmap(nA,'CellLabelColor','none');
h.FontSize = 17;
h.Colormap = cool(4);
h.GridVisible = 'off';
set(gcf, 'Position',  [10 10 700 600])
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
exportgraphics(gca, 'figures/normsaylr3.pdf')

% gca = plot_hmat_nprec(aphA);
% exportgraphics(h, 'figures/pnormLeGres.pdf')

data = compute_hmat_norm(aphA);
pA = compute_hmat_prec(aphA);
% h = heatmap(nA);
ColorMap = cool(4);

[m, n] = size(data);
fig = figure('Renderer', 'painters', 'Position', [10 10 750 600]);
ax = axes(fig);
h = imagesc(ax, data);
set(ax,'XTick',1:m,'YTick',1:n);
% title('imagesc')
ax.TickLength(1) = 0;
% Create heatmap's colormap
cmap = ColorMap; % [linspace(.9,0,n)', linspace(.9447,.447,n)', linspace(.9741,.741,n)'];
h2 = colormap(ax, cmap); 
h3 = colorbar(ax);

hold on
% Set grid lines
% arrayfun(@(x)xline(ax,x,'k-','Alpha',1),0.5:1:(m+.5));
% arrayfun(@(y)yline(ax,y,'k-','Alpha',1),0.5:1:(n+.5));
[linesF, columnsF, valuesF] = find(pA);
% th = text(linesF, columnsF, string(valuesF), ...
%    'VerticalAlignment', 'middle','HorizontalAlignment','Center');
fontsize(9,"points");
exportgraphics(gca, 'figures/pnormsaylr3.pdf');