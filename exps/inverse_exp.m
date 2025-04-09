% Add mhodlr to path
addpath("../mhodlr/");
% Set random seed for reproducibility
rng(0);

n = 2000; % Matrix size
depths = [5]; % Depths to test
vareps_values = [1e-8, 1e-12]; % Approximation tolerances
precisions = {'s', 'b'}; % Precision types (fp32, half, b)
precision_labels = {'fp32', 'b'}; % For plot labels
min_block_size = 10; % Minimum block size
method = 'svd'; % Compression method
matrix_names = {'mat-1', 'mat-2', 'mat-3', 'mat-4'}; % Kernel matrix labels

% Generate point sets
% s1: 1D uniform grid [0, 1]
s1 = linspace(0, 1, n);

% s2: 2D uniform grid [-1, 1] x [-1, 1], approximately 2000 points
grid_size = ceil(sqrt(n)); % e.g., 45 for n=2000
[x, y] = meshgrid(linspace(-1, 1, grid_size), linspace(-1, 1, grid_size));
s2 = [x(:), y(:)]; % Flatten to list of points
s2 = s2(1:n, :); % Trim to exactly 2000 points

% Function to generate kernel matrices
generate_kernel = @(type, points, h) generate_kernel_matrix(type, points, h);

% Initialize storage for all matrices
errors = cell(length(matrix_names), length(depths)); % Matrix x Depth
for m = 1:length(matrix_names)
    for d = 1:length(depths)
        errors{m, d} = zeros(length(vareps_values), length(precisions));
    end
end

% Simulate for each kernel matrix
for m = 1:length(matrix_names)
    fprintf('Processing %s\n', matrix_names{m});
    
    % Generate kernel matrix
    if m == 1 % mat-1: Kernel (i) on s1
        A = generate_kernel(1, s1, []);
    elseif m == 2 % mat-2: Kernel (ii) on s2
        A = generate_kernel(2, s2, []);
    elseif m == 3 % mat-3: Kernel (iii) on s2, h=1
        A = generate_kernel(3, s2, 1);
    elseif m == 4 % mat-4: Kernel (iii) on s2, h=20
        A = generate_kernel(3, s2, 20);
    end
    
    % Compute true inverse once (in double precision)
    A_inv = inv(A); % MATLAB's dense inverse
    norm_A_inv = norm(A_inv, 'fro'); % For relative error normalization
    
    % Loop over depths
    for d = 1:length(depths)
        depth = depths(d);
        fprintf('  Depth = %d\n', depth);
        
        % Loop over vareps
        for v = 1:length(vareps_values)
            vareps = vareps_values(v);
            
            % Construct HODLR matrix
            hA = hodlr(A, depth, min_block_size, method, vareps, 99999);
            
            % Loop over precisions
            for p = 1:length(precisions)
                prec = precisions{p};
                u = precision(prec);
                set_prec(u);
                
                % Compute HODLR inverse
                iA = minverse(hA);
                
                % Convert to dense and compute relative error
                iA_dense = iA.todense(); % Assuming todense() recovers dense matrix
                error = norm(iA_dense - A_inv, 'fro') / norm_A_inv;
                errors{m, d}(v, p) = error;
                fprintf('    vareps %e, precision %s: relative error = %e\n', ...
                    vareps, precision_labels{p}, error);
            end
        end
    end
end

% Plotting and saving figures
for m = 1:length(matrix_names)
    for d = 1:length(depths)
        figure;
        bar(errors{m, d});
        set(gca, 'XTickLabel', cellstr(num2str(vareps_values', '%.0e')));
        set(gca, 'YScale', 'log'); % Logarithmic y-axis
        xlabel('vareps');
        ylabel('Relative Error (Frobenius Norm)');
        title(sprintf('%s, Depth = %d', matrix_names{m}, depths(d)));
        legend(precision_labels, 'Location', 'best');
        grid on;
        
        % Save as JPG
        filename = sprintf('%s_depth%d.jpg', matrix_names{m}, depths(d));
        saveas(gcf, filename, 'jpg');
        close(gcf); % Close figure to free memory
    end
end

% Helper function to generate kernel matrices
function A = generate_kernel_matrix(type, points, h)
    n = size(points, 1);
    A = zeros(n, n);
    
    if type == 1 % Kernel (i): 1/(x - y)
        x = points(:); % 1D points
        for i = 1:n
            for j = 1:n
                if i == j
                    A(i, j) = 1;
                else
                    A(i, j) = 1 / (x(i) - x(j));
                end
            end
        end
    elseif type == 2 % Kernel (ii): log ||x - y||_2
        for i = 1:n
            for j = 1:n
                if i == j
                    A(i, j) = 0;
                else
                    diff = points(i, :) - points(j, :);
                    A(i, j) = log(norm(diff, 2));
                end
            end
        end
    elseif type == 3 % Kernel (iii): exp(-||x - y||_2^2 / (2h^2))
        for i = 1:n
            for j = 1:n
                diff = points(i, :) - points(j, :);
                A(i, j) = exp(-norm(diff, 2)^2 / (2 * h^2));
            end
        end
    end
end