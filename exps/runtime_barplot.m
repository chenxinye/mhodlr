
hodlr_data = readtable('hmtoolbox_runtimes.csv');
mhodlr_data = readtable('mhodlr_runtimes.csv');

matrix_sizes = hodlr_data.MatrixSize; % Assuming both files have same sizes
hodlr_construction = hodlr_data.ConstructionTime;
mhodlr_construction = mhodlr_data.ConstructionTime;
hodlr_matvec = hodlr_data.MatVecTime;
mhodlr_matvec = mhodlr_data.MatVecTime;

% Plot 1: Construction Runtime Comparison
figure('Name', 'HODLR vs mHODLR Construction Runtime');
bar_data_construction = [hodlr_construction, mhodlr_construction];
bar(matrix_sizes, bar_data_construction, 'grouped');
xlabel('Matrix Size');
ylabel('Construction Time (seconds)');
title('Construction Runtime Comparison: HODLR vs mHODLR');
legend('hm-toolbox', 'mHODLR', 'Location', 'northwest');
grid on;
set(gcf, 'Position', [100, 100, 600, 400]);
print('-djpeg', 'construction_comparison.jpg'); 
close; 


% Plot 2: Matrix-Vector Runtime Comparison
figure('Name', 'HODLR vs mHODLR Matrix-Vector Runtime');
bar_data_matvec = [hodlr_matvec, mhodlr_matvec];
bar(matrix_sizes, bar_data_matvec, 'grouped');
xlabel('Matrix Size');
ylabel('Matrix-Vector Time (seconds)');
title('Matrix-Vector Runtime Comparison: HODLR vs mHODLR');
legend('hm-toolbox', 'mHODLR', 'Location', 'northwest');
grid on;
set(gcf, 'Position', [100, 100, 600, 400]); 
print('-djpeg', 'matvec_comparison.jpg');
close; 

% set(gcf, 'Position', [100, 100, 600, 400]); % Adjust size of the second figure