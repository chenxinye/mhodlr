function [VA, h] = plot_hmat_norm(obj, VA, varargin)
    
    if nargin == 2
        level = 1;
        VA = repmat(obj.normOrder(1), [2^obj.bottom_level, 2^obj.bottom_level]);
    else
        level = varargin{1};
    end

    [rowSize, colSize] = size(VA);

    if level > obj.bottom_level | level > obj.max_level
        VA(:, :) = obj.normOrder(obj.level+1);
        return;
    else
        obj.level = level;
        level = level + 1;
        rowSplit = ceil(rowSize / 2);
        colSplit = ceil(colSize / 2);
        VA(1:rowSplit, colSplit+1:end) = obj.normOrder(obj.level+1);
        VA(rowSplit+1:end, 1:colSplit) = obj.normOrder(obj.level+1);
        VA(1:rowSplit, 1:colSplit) = plot_hmat_norm(obj, VA(1:rowSplit, 1:colSplit), level);
        VA(rowSplit+1:end, colSplit+1:end) = plot_hmat_norm(obj, VA(rowSplit+1:end, colSplit+1:end), level);
        
    end

    h = heatmap(VA);
    % h.GridVisible = 'off';

end
