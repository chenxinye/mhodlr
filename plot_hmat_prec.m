function [VA, h] =  plot_hmat_prec(obj, VA, precs, level)
    
    if level == 1
        VA = zeros(2^obj.bottom_level, 2^obj.bottom_level);
    end

    [rowSize, colSize] = size(VA);

    if level > obj.bottom_level | level > (obj.max_level-1)
        VA(:, :) = precs(obj.level);
        return;
    else
        obj.level = level;

        level = level + 1;
        rowSplit = ceil(rowSize / 2);
        colSplit = ceil(colSize / 2);

        VA(1:rowSplit, colSplit+1:end) = precs(obj.level);
        VA(rowSplit+1:end, 1:colSplit) = precs(obj.level);

        VA(1:rowSplit, 1:colSplit) = plot_hmat_prec(obj, VA(1:rowSplit, 1:colSplit), precs, level);
        VA(rowSplit+1:end, colSplit+1:end) = plot_hmat_prec(obj, VA(rowSplit+1:end, colSplit+1:end), precs, level);
        
    end

    h = heatmap(VA);
    % h.GridVisible = 'off';

end
