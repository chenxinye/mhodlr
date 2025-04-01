function visualize_hodlr(H, level, colors)
    % VISUALIZE_HODLR - Visualize HODLR matrix partitions with user-specified colors
    %   H: hodlr object from mhodlr library
    %   colors: cell array of color specs (e.g., {'r', 'g', 'b', 'y'}), one per level or block
    %   level: current recursion level (optional, default = 1)
    
    if nargin < 2 || isempty(colors)
        % Default colors if none provided
        colors = {'r', 'g', 'b', 'y', 'c', 'm', 'k'};
    end
    if nargin < 3
        level = 1;  % Start at top level
    end
    
    % Ensure figure is ready
    figure;
    hold on;
    axis equal;
    axis off;
    title(['HODLR Matrix Visualization (Level ', num2str(level), ')']);
    
    % Total matrix dimensions
    m = H.shape(1);
    n = H.shape(2);
    
    % Recursively plot the HODLR structure
    plot_hodlr_block(H, 0, 0, m, n, colors, level, 1);
    
    % Set plot limits
    xlim([0 n]);
    ylim([0 m]);
    set(gca, 'YDir', 'reverse');  % Matrix convention: (0,0) at top-left
    hold off;
end

function plot_hodlr_block(H, x0, y0, m, n, colors, level, max_level)
    % PLOT_HODLR_BLOCK - Recursively plot HODLR blocks
    %   H: hodlr object
    %   x0, y0: bottom-left corner of current block
    %   m, n: rows and columns of current block
    %   colors: user-specified colors
    %   level: current level
    %   max_level: maximum depth from top (tracked for color assignment)
    
    % Select color based on level (cycle through colors if needed)
    color_idx = mod(level - 1, length(colors)) + 1;
    block_color = colors{color_idx};
    
    if ~isempty(H.D)
        % Leaf node: plot dense block
        rectangle('Position', [x0, y0, n, m], ...
                  'EdgeColor', block_color, 'LineWidth', 2, ...
                  'FaceColor', [block_color, 0.2]);  % Slight transparency
        text(x0 + n/2, y0 + m/2, 'D', 'HorizontalAlignment', 'center', ...
             'Color', 'k', 'FontWeight', 'bold');
    else
        % Non-leaf node: compute block sizes
        sv1 = size(H.V1, 2);  % Columns of A11
        su1 = size(H.U1, 1);  % Rows of A11
        m2 = m - su1;         % Rows of A22
        n2 = n - sv1;         % Columns of A22
        
        % Plot A11 (top-left diagonal block)
        plot_hodlr_block(H.A11, x0, y0 + m2, su1, sv1, colors, level + 1, max_level);
        
        % Plot A22 (bottom-right diagonal block)
        plot_hodlr_block(H.A22, x0 + sv1, y0, m2, n2, colors, level + 1, max_level);
        
        % Plot U1 V2 (top-right off-diagonal block)
        rectangle('Position', [x0 + sv1, y0 + m2, n2, su1], ...
                  'EdgeColor', block_color, 'LineWidth', 1.5, ...
                  'FaceColor', [block_color, 0.1]);
        text(x0 + sv1 + n2/2, y0 + m2 + su1/2, 'U1 V2', ...
             'HorizontalAlignment', 'center', 'Color', 'k');
        
        % Plot U2 V1 (bottom-left off-diagonal block)
        rectangle('Position', [x0, y0, sv1, m2], ...
                  'EdgeColor', block_color, 'LineWidth', 1.5, ...
                  'FaceColor', [block_color, 0.1]);
        text(x0 + sv1/2, y0 + m2/2, 'U2 V1', ...
             'HorizontalAlignment', 'center', 'Color', 'k');
    end
end

% Helper function to mimic is_hodlr_class (for completeness)
function is_hodlr = is_hodlr_class(obj)
    switch class(obj)
        case {'hodlr', 'amphodlr', 'mphodlr'}
            is_hodlr = true;
        otherwise
            is_hodlr = false;
    end
end