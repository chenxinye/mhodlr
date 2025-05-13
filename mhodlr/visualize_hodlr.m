function visualize_hodlr(H, colors, output_filename)
%{
    VISUALIZE_HODLR - Visualize HODLR matrix partitions with dynamic colors and save figure
        Parameters
        --------------------
        H - hodlr, mphodlr, and amphodlr
            hodlr object from mhodlr library

        colors - cell array
            optional cell array of color specs (e.g., {'r', 'g', 'b'}), defaults to auto-generated colors based on bottom_level

        output_filename - string
            optional filename for saving the figure (e.g., 'hodlr_plot.png'). If not provided, defaults to 'hodlr_visualization.png'.
%}
   
    if ~is_hodlr_class(H)
        error('Input must be a hodlr, amphodlr, or mphodlr object.');
    end
    
    % Get bottom level (depth) from the hodlr object
    bottom_level = H.bottom_level;  % Maximum depth (e.g., 3 for levels 0,1,2,3)
    num_levels = bottom_level + 1;  % Number of levels (0 to bottom_level)
    
    % Generate default colors if not provided
    if nargin < 2 || isempty(colors)
        color_matrix = lines(num_levels + 1);  % +1 for distinct dense block color
        colors = num2cell(color_matrix, 2);    % Convert rows to cell array
        for i = 1:num_levels + 1
            colors{i} = color_matrix(i, :);
        end
    elseif length(colors) < num_levels + 1
        warning('Number of colors (%d) is less than required (%d). Colors will cycle.', ...
                length(colors), num_levels + 1);
    end
    
    % Set default output filename if not provided
    if nargin < 3 || isempty(output_filename)
        output_filename = 'hodlr_visualization.png';
    end
    
    % Split colors: first num_levels for levels, last for dense blocks
    level_colors = colors(1:num_levels);  % Colors for each level (0 to bottom_level)
    dense_color = colors{end};            % Distinct color for dense blocks (D)
    
    % Ensure figure is ready
    fig = figure;
    hold on;
    axis equal;
    axis off;
    title('HODLR Structure Visualization');
    
    % Total matrix dimensions
    m = H.shape(1);
    n = H.shape(2);
    
    % Recursively plot the HODLR structure
    plot_hodlr_block(H, 0, 0, m, n, level_colors, dense_color);
    
    % Set plot limits
    xlim([0 n]);
    ylim([0 m]);
    set(gca, 'YDir', 'reverse');  % Matrix convention: (0,0) at top-left
    hold off;
    
    % Save the figure to the current working directory
    saveas(fig, fullfile(pwd, output_filename));
end

function plot_hodlr_block(H, x0, y0, m, n, level_colors, dense_color)
    % PLOT_HODLR_BLOCK - Recursively plot HODLR blocks with corrected left-right orientation
    %   H: hodlr object
    %   x0, y0: bottom-left corner of current block
    %   m, n: rows and columns of current block
    %   level_colors: cell array of colors for each level
    %   dense_color: distinct color for dense blocks (D)
    
    % Get current level from the hodlr object
    current_level = H.level;  
    % Assumes this field exists in mhodlr
    
    % Select color for this level (cycle if needed)
    color_idx = mod(current_level, length(level_colors)) + 1;
    block_color = level_colors{color_idx};
    
    if ~isempty(H.D)
        % Leaf node: plot dense block with actual dimensions from H.D
        d_m = size(H.D, 1);  % Rows of D
        d_n = size(H.D, 2);  % Columns of D
        rectangle('Position', [x0, y0, d_n, d_m], ...
                  'EdgeColor', dense_color, 'LineWidth', 2, ...
                  'FaceColor', [dense_color, 0.2]);  % Slight transparency
        text(x0 + d_n/2, y0 + d_m/2, 'D', 'HorizontalAlignment', 'center', ...
             'Color', 'k', 'FontWeight', 'bold');
    else
        % Non-leaf node: compute block sizes
        sv1 = size(H.V1, 2);  % Columns of A11
        su1 = size(H.U1, 1);  % Rows of A11
        m2 = m - su1;         % Rows of A22
        n2 = n - sv1;         % Columns of A22
        
        % Plot A11 (top-left diagonal block, adjust x0 to right for reversal)
        plot_hodlr_block(H.A11, x0 + n2, y0 + m2, su1, sv1, level_colors, dense_color);
        
        % Plot A22 (bottom-right diagonal block, adjust x0 to left)
        plot_hodlr_block(H.A22, x0, y0, m2, n2, level_colors, dense_color);
        
        % Plot U1 V2 (top-right off-diagonal block, move to top-left)
        rectangle('Position', [x0, y0 + m2, n2, su1], ...
                  'EdgeColor', block_color, 'LineWidth', 1.5, ...
                  'FaceColor', [block_color, 0.1]);
        text(x0 + n2/2, y0 + m2 + su1/2, 'U1 V2', ...
             'HorizontalAlignment', 'center', 'Color', 'k');
        
        % Plot U2 V1 (bottom-left off-diagonal block, move to bottom-right)
        rectangle('Position', [x0 + n2, y0, sv1, m2], ...
                  'EdgeColor', block_color, 'LineWidth', 1.5, ...
                  'FaceColor', [block_color, 0.1]);
        text(x0 + n2 + sv1/2, y0 + m2/2, 'U2 V1', ...
             'HorizontalAlignment', 'center', 'Color', 'k');
    end
end

function is_hodlr = is_hodlr_class(obj)
    switch class(obj)
        case {'hodlr', 'amphodlr', 'mphodlr'}
            is_hodlr = true;
        otherwise
            is_hodlr = false;
    end
end
