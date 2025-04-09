addpath("../mhodlr/");

rng(0);

n = 200; % Matrix size
depths = [3, 8]; % Depths to test
vareps_values = [1e-4, 1e-8, 1e-12]; % Approximation tolerances
precisions = {'s', 'h', 'b'}; % Precision types (fp32, half, b)
precision_labels = {'fp32', 'half', 'b'}; % For CSV column headers
min_block_size = 10; % Minimum block size
method = 'svd'; % Compression method
matrix_names = {'mat-1', 'mat-2', 'mat-3', 'mat-4'}; % Kernel matrix labels
lambda = 1; % Initial diagonal shift for SPD enforcement

% Generate point sets
s1 = linspace(0, 1, n); % 1D uniform grid [0, 1]
grid_size = ceil(sqrt(n)); % e.g., 45 for n=2000
[x, y] = meshgrid(linspace(-1, 1, grid_size), linspace(-1, 1, grid_size));
s2 = [x(:), y(:)]; % 2D uniform grid [-1, 1] x [-1, 1]
s2 = s2(1:n, :); % Trim to 2000 points

% Function to generate kernel matrices
generate_kernel = @(type, points, h) generate_kernel_matrix(type, points, h);

% Initialize storage
errors = cell(length(matrix_names), length(depths));
for m = 1:length(matrix_names)
    for d = 1:length(depths)
        errors{m, d} = zeros(length(vareps_values), length(precisions));
    end
end

% Simulate for each kernel matrix
for m = 1:length(matrix_names)
    fprintf('Processing %s\n', matrix_names{m});
    
    % Generate kernel matrix and random vector
    try
        if m == 1 % mat-1: Modified Kernel (i) on s1
            A = generate_kernel(1, s1, []);
            fprintf('  Note: mat-1 modified to be SPD.\n');
        elseif m == 2 % mat-2: Kernel (ii) on s2
            A = generate_kernel(2, s2, []);
            fprintf('  Note: mat-2 shifted to be SPD.\n');
        elseif m == 3 % mat-3: Kernel (iii) on s2, h=1
            A = generate_kernel(3, s2, 1);
            fprintf('  Note: mat-3 is naturally SPD.\n');
        elseif m == 4 % mat-4: Kernel (iii) on s2, h=20
            A = generate_kernel(3, s2, 20);
            fprintf('  Note: mat-4 is naturally SPD.\n');
        end
        
        if any(isnan(A(:))) || any(isinf(A(:)))
            error('NaN or Inf detected in kernel matrix');
        end
        
        % Enforce SPD with diagonal shift
        shift = lambda;
        A = A + shift * eye(n); % Add initial shift
        spd_verified = false;
        while ~spd_verified
            try
                chol(A); % Test if SPD
                spd_verified = true;
            catch
                fprintf('  %s not SPD with shift %f, increasing shift\n', matrix_names{m}, shift);
                shift = shift * 10; % Increase shift if not SPD
                A = A - (shift/10) * eye(n) + shift * eye(n); % Update with larger shift
            end
        end
        fprintf('  %s made SPD with shift %f\n', matrix_names{m}, shift);
        
        % Generate random vector and true product
        x = rand(n, 1);
        Ax = A * x;
        norm_Ax = norm(Ax, 'fro');
        if isnan(norm_Ax) || isinf(norm_Ax)
            error('Invalid norm of A*x');
        end
    catch e
        fprintf('Error generating %s: %s\n', matrix_names{m}, e.message);
        for d = 1:length(depths)
            errors{m, d}(:) = -1;
        end
        continue;
    end
    
    % Loop over depths
    for d = 1:length(depths)
        depth = depths(d);
        fprintf('  Depth = %d\n', depth);
        
        % Loop over vareps
        for v = 1:length(vareps_values)
            vareps = vareps_values(v);
            
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
                    u = precision(prec);
                    set_prec(u);
                    b = hdot(hA, x, 'dense'); % HODLR matrix-vector product
                    
                    % Compute relative error
                    if any(isnan(b(:))) || any(isinf(b(:)))
                        error('NaN or Inf in hA*x');
                    end
                    error = norm(b - Ax, 'fro') / norm_Ax;
                    if isnan(error) || isinf(error)
                        error('Invalid error value');
                    end
                catch e
                    fprintf('    vareps %e, precision %s: Error (%s)\n', ...
                        vareps, precision_labels{p}, e.message);
                    error = -1;
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
        T = array2table(errors{m, d}, 'VariableNames', precision_labels);
        T.vareps = vareps_values';
        T = movevars(T, 'vareps', 'Before', 1);
        filename = sprintf('mv_%s_depth%d_matvec.csv', matrix_names{m}, depths(d));
        writetable(T, filename);
        fprintf('Saved results to %s\n', filename);
    end
end

% Helper function to generate kernel matrices (modified for SPD)
function A = generate_kernel_matrix(type, points, h)
    n = size(points, 1);
    A = zeros(n, n);
    
    if type == 1 % Modified Kernel (i): 1/|x - y|
        x = points(:);
        for i = 1:n
            for j = 1:n
                if i == j
                    A(i, j) = 1; % Keep diagonal as 1
                else
                    A(i, j) = 1 / abs(x(i) - x(j)); % Symmetric
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