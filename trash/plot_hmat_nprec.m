function gca = plot_hmat_nprec(obj, varargin)
    if nargin == 1
        colorMap = spring;
    else if nargin == 2
        colorMap = varargin{1};
    else
        colorMap = varargin{1};
        colorMap = colorMap(varargin{2});
    end

    data = compute_hmat_norm(obj);
    pA = compute_hmat_prec(obj);
    % h = heatmap(nA);
    % h.Colormap = parula;
    
    [m, n] = size(data)
    fig = figure('Renderer', 'painters', 'Position', [10 10 750 600])
    ax = axes(fig);
    h = imagesc(ax, data);
    set(ax,'XTick',1:m,'YTick',1:n)
    % title('imagesc')
    ax.TickLength(1) = 0;
    % Create heatmap's colormap
    cmap = colorMap; % [linspace(.9,0,n)', linspace(.9447,.447,n)', linspace(.9741,.741,n)'];
    h2 = colormap(ax, cmap); 
    h3 = colorbar(ax)
    hold on
    % Set grid lines
    arrayfun(@(x)xline(ax,x,'k-','Alpha',1),0.5:1:(m+.5))
    arrayfun(@(y)yline(ax,y,'k-','Alpha',1),0.5:1:(n+.5))
    [linesF, columnsF, valuesF] = find(pA);
    th = text(linesF, columnsF, string(valuesF), ...
        'VerticalAlignment', 'middle','HorizontalAlignment','Center');
    
    % h.GridVisible = 'off';
end
