function [VA, h] =  plot_hmat_prec(obj, VA, varargin)
    
    if nargin == 2
        level = 1;
        VA = zeros(2^obj.bottom_level, 2^obj.bottom_level);
    else
        level = varargin{1};
    end

    [rowSize, colSize] = size(VA);

    if level > obj.bottom_level | level > obj.max_level
        VA(:, :) = obj.precIndex(obj.level);
        return;
    else
        obj.level = level;

        level = level + 1;
        rowSplit = ceil(rowSize / 2);
        colSplit = ceil(colSize / 2);

        VA(1:rowSplit, colSplit+1:end) = obj.precIndex(obj.level);
        VA(rowSplit+1:end, 1:colSplit) = obj.precIndex(obj.level);

        VA(1:rowSplit, 1:colSplit) = plot_hmat_prec(obj, VA(1:rowSplit, 1:colSplit), level);
        VA(rowSplit+1:end, colSplit+1:end) = plot_hmat_prec(obj, VA(rowSplit+1:end, colSplit+1:end), level);
        
    end

    h = heatmap(VA);
    h.Colormap = parula;
    h.ColorbarVisible = 'off';
    % h.GridVisible = 'off';

end
