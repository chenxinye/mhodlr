% Parameters
matrix_names = {'mat-1', 'mat-2', 'mat-3', 'mat-4'};
depths = [3, 8];
vareps_values = [1e-4, 1e-8, 1e-12];
precision_labels = {'fp32', 'tf32', 'half', 'bf16'}; % For CSV column headers

% Colors for each vareps
colors = {'b', 'g', 'r'}; % Blue, Green, Red for 1e-4, 1e-8, 1e-12

% Loop over matrices and depths
for m = 1:length(matrix_names)
    for d = 1:length(depths)
        % Read CSV file
        filename = sprintf('lu_%s_depth%d.csv', matrix_names{m}, depths(d));
        if ~exist(filename, 'file')
            fprintf('File %s not found, skipping.\n', filename);
            continue;
        end
        data = readtable(filename);
        
        % Extract error data (excluding vareps column)
        errors = table2array(data(:, precision_labels));
        
        % Create figure
        figure('Position', [100, 100, 800, 600]);
        hold on;
        
        % Bar plot setup
        bar_width = 0.2;
        x = 1:length(precision_labels); % x-axis positions for precisions
        
        % Plot bars for each vareps
        for v = 1:length(vareps_values)
            bar_pos = x + (v-2) * bar_width; % Offset for each vareps
            bar(bar_pos, errors(v, :), bar_width, 'FaceColor', colors{v}, ...
                'DisplayName', sprintf('vareps = %e', vareps_values(v)));
        end
        
        % Customize plot
        set(gca, 'XTick', x, 'XTickLabel', precision_labels, 'FontSize', 12);
        xlabel('Precision Type', 'FontSize', 14);
        ylabel('Relative Error', 'FontSize', 14);
        title(sprintf('%s, Depth %d (LU)', matrix_names{m}, depths(d)), 'FontSize', 16);
        legend('Location', 'best', 'FontSize', 12);
        
        % Set y-axis to log scale for better visibility of small errors
        set(gca, 'YScale', 'log');
        ylim([1e-15, 1]); % Adjust based on your error range
        
        % Replace -1 with a small value for visibility (if any errors are -1)
        errors(errors == -1) = 1e-15; % Below ylim for failed runs
        
        hold off;
        
        % Save plot
        plot_filename = sprintf('lu_%s_depth%d_barplot.png', matrix_names{m}, depths(d));
        saveas(gcf, plot_filename);
        fprintf('Saved bar plot to %s\n', plot_filename);
        
        % Close figure to free memory
        close(gcf);
    end
end