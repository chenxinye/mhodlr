% Add mhodlr to path
addpath("../mhodlr/");

% Set random seed for reproducibility
rng(2025);

% Parameters
n = 100; % Matrix size
depths = [3, 8]; % Depths to test
vareps_values = [1e-4, 1e-8, 1e-12]; % Approximation tolerances
precisions = {'s', 'h', 't', 'b'}; % Precision types (fp32, half)
precision_labels = {'fp32', 'half', 'tf32', 'b'}; % For CSV column headers
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
    try
        if m == 1 % mat-1: Kernel (i) on s1
            A = generate_kernel(1, s1, []);
        elseif m == 2 % mat-2: Kernel (ii) on s2
            A = generate_kernel(2, s2, []);
        elseif m == 3 % mat-3: Kernel (iii) on s2, h=1
            A = generate_kernel(3, s2, 1);
        elseif m == 4 % mat-4: Kernel (iii) on s2, h=20
            A = generate_kernel(3, s2, 20);
        end
        
        % Check for NaN in A
        if any(isnan(A(:)))
            error('NaN detected in kernel matrix');
        end
        
        % Compute norm of A for relative error
        norm_A = norm(A, 'fro');
        if isnan(norm_A) || isinf(norm_A)
            error('Invalid norm of A');
        end
    catch e
        fprintf('Error generating %s: %s\n', matrix_names{m}, e.message);
        for d = 1:length(depths)
            errors{m, d}(:) = -1; % Set all errors to -1 if matrix generation fails
        end
        continue; % Skip to next matrix
    end
    
    % Loop over depths
    for d = 1:length(depths)
        depth = depths(d);
        fprintf('  Depth = %d\n', depth);
        
        % Loop over vareps
        for v = 1:length(vareps_values)
            vareps = vareps_values(v);
            
            % Construct HODLR matrix with exception handling
            try
                hA = hodlr(A, depth, min_block_size, method, vareps);
            catch e
                fprintf('    vareps %e: HODLR construction failed (%s)\n', vareps, e.message);
                errors{m, d}(v, :) = -1;
                continue;
            end
            
            % Loop over precisions
            for p = 1:length(precisions)
                prec = precisions{p};
                try
                    if strcmp(prec, 's') % fp32: use full precision hlu
                        [L, U] = hlu(hA);
                    else % half: use mixed-precision mhlu
                        u = precision(prec);
                        set_prec(u);
                        [L, U] = mhlu(hA);
                    end
                    
                    % Compute relative error with checks
                    LU_product = hdot(L, U, 'dense');
                    if any(isnan(LU_product(:))) || any(isinf(LU_product(:)))
                        error('NaN or Inf in LU product');
                    end
                    error = norm(LU_product - A, 'fro') / norm_A;
                    if isnan(error) || isinf(error)
                        error('Invalid error value');
                    end
                catch e
                    fprintf('    vareps %e, precision %s: Error (%s)\n', ...
                        vareps, precision_labels{p}, e.message);
                    error = -1; % Set to -1 on failure
                end
                errors{m, d}(v, p) = error;
                fprintf('    vareps %e, precision %s: relative error = %e\n', ...
                    vareps, precision_labels{p}, error);
            end
        end
    end
end

% Save results to CSV files
for m = 1:length(matrix_names)
    for d = 1:length(depths)
        % Create table with vareps as rows and precisions as columns
        T = array2table(errors{m, d}, 'VariableNames', precision_labels);
        T.vareps = vareps_values'; % Add vareps as a column
        T = movevars(T, 'vareps', 'Before', 1); % Move vareps to first column
        
        % Save to CSV
        filename = sprintf('%s_depth%d.csv', matrix_names{m}, depths(d));
        writetable(T, filename);
        fprintf('Saved results to %s\n', filename);
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