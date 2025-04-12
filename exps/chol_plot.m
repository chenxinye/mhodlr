% Parameters
matrix_names = {'mat-1', 'mat-2', 'mat-3', 'mat-4'};
depths = [3, 8];
vareps_values = [1e-2, 1e-5, 1e-8, 1e-11];
precision_labels = {'fp32', 'fp16', 'bf16'};
operation = 'cholesky';

% Colors for each vareps (RGB for consistency)
colors = [
    0, 0, 1;    % Blue for 1e-2
    0, 0.5, 0;  % Green for 1e-5
    1, 0, 0;    % Red for 1e-8
    1, 1, 0     % Yellow for 1e-11
];

% Loop over matrices
for m = 1:length(matrix_names)
    % Create figure
    figure('Position', [100, 100, 1200, 500]);
    
    % Store handles for legend
    legend_handles = [];
    legend_labels = {};
    
    for d = 1:length(depths)
        % Read CSV file
        filename = sprintf('chol_%s_depth%d.csv', matrix_names{m}, depths(d));
        if ~exist(filename, 'file')
            fprintf('File %s not found, skipping depth %d for %s.\n', filename, depths(d), matrix_names{m});
            continue;
        end
        data = readtable(filename);
        
        % Extract error data (excluding vareps column)
        errors = table2array(data(:, precision_labels));
        
        % Create axes
        axes('Position', [0.05 + (d-1)*0.45, 0.2, 0.4, 0.7]);
        hold on;
        
        % Bar plot setup
        bar_width = 0.2;
        x = 1:length(precision_labels); % x-axis positions for precisions
        
        % Plot bars for each vareps
        for v = 1:length(vareps_values)
            bar_pos = x + (v-2) * bar_width; % Offset for each vareps
            h = bar(bar_pos, errors(v, :), bar_width, 'FaceColor', colors(v, :));
            % Store handle and label for first depth only
            if d == 1
                legend_handles(end+1) = h;
                % Use scientific notation: vareps = 10^-X
                exponent = log10(vareps_values(v));
                legend_labels{end+1} = sprintf('vareps = 10^{%d}', exponent);
            end
        end
        
        % Customize axes
        set(gca, 'XTick', x, 'XTickLabel', precision_labels, 'FontSize', 14);
        xlabel('Precision Type', 'FontSize', 12);
        ylabel('Relative Error', 'FontSize', 12);
        title(sprintf('Depth %d', depths(d)), 'FontSize', 14);
        
        % Set y-axis to log scale
        set(gca, 'YScale', 'log');
        ylim([1e-15, 1]); % Adjust based on your error range
        
        % Replace -1 with a small value for visibility
        errors(errors == -1) = 1e-15; % Below ylim for failed runs
        
        hold off;
    end
    
    % Add super title
    sgtitle(sprintf('%s (Cholesky)', matrix_names{m}), 'FontSize', 16);
    
    % Add shared legend below subplots
    legend(legend_handles, legend_labels, 'Orientation', 'horizontal', ...
           'Location', 'southoutside', 'FontSize', 14, ...
           'Position', [0.25, 0.05, 0.5, 0.05]);
    
    % Save plot
    plot_filename = sprintf('%s_%s_barplot.png', matrix_names{m}, operation);
    saveas(gcf, plot_filename);
    fprintf('Saved bar plot to %s\n', plot_filename);
    
    % Close figure to free memory
    close(gcf);
end