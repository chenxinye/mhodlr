
addpath('../mhodlr');

matrix_sizes = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000];
min_block_size = 200;
vareps = 1e-6;
num_trials = 5;

construction_times = zeros(length(matrix_sizes), 1);
matvec_times = zeros(length(matrix_sizes), 1);

% Loop over matrix sizes
for i = 1:length(matrix_sizes)
    n = matrix_sizes(i);
    fprintf('Testing matrix size: %d\n', n);
    
    % Generate a sample kernel matrix (e.g., 1/|x_i - x_j|)
    x = linspace(0, 1, n)';
    A = rand(n); %zeros(n, n); % Preallocate dense matrix
    % for row = 1:n
    %     for col = 1:n
    %         A(row, col) = 1 / (abs(x(row) - x(col)) + 1e-10); % Evaluate function
    %     end
    % end
    
    H = hodlr(A, 999, min_block_size, "svd", vareps);
    v = randn(n, 1);
    w = hdot(H, v, 'dense'); % Warm-up mat-vec
    
    % Measure construction time (average over trials after warm-up)
    constr_time = 0;
    for t = 1:(num_trials - 1) % One less than total since warm-up is discarded
        tic;
        H = hodlr(A, 999, min_block_size, "svd", vareps);
        constr_time = constr_time + toc;
    end
    
    construction_times(i) = constr_time / (num_trials - 1);
    
    % Generate a random vector for matrix-vector product
    v = randn(n, 1);
    
    % Measure matrix-vector product time (average after warm-up)
    matvec_time = 0;
    for t = 1:(num_trials - 1)
        tic;
        w = hdot(H, v, 'dense');
        matvec_time = matvec_time + toc;
    end
    matvec_times(i) = matvec_time / (num_trials - 1);
end

results = table(matrix_sizes', construction_times, matvec_times, ...
    'VariableNames', {'MatrixSize', 'ConstructionTime', 'MatVecTime'});

filename = 'mhodlr_runtimes.csv';
writetable(results, filename);
fprintf('Results saved to %s\n', filename);

disp(results);