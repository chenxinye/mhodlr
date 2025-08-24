function [R, S] = two_sided_diagonal_scaling_sym(A, tol, max_iter)
    % Symmetry-preserving two-sided diagonal scaling for a symmetric matrix A.
    %
    % Parameters:
    %   A (matrix): Input symmetric matrix
    %   tol (double): Convergence tolerance (default: 1e-6)
    %   max_iter (int): Maximum number of iterations (default: 100)
    %
    % Returns:
    %   R (vector): Row scaling vector
    %   S (vector): Column scaling vector
    if nargin < 2
        tol = 1e-6;
    end
    if nargin < 3
        max_iter = 2;
    end
    
    [m, n] = size(A);
    
    % Initialize scaling vectors
    R = ones(m, 1);
    S = ones(n, 1);
    
    for iter = 1:max_iter
        % Compute row and column sums efficiently
        row_sums = sqrt(sum(abs(A), 2));
        col_sums = sqrt(sum(abs(A), 1));
        
        % Avoid division by zero
        row_sums(row_sums == 0) = 1;
        col_sums(col_sums == 0) = 1;
        
        % Compute scaling factors
        R_new = 1 ./ row_sums;
        S_new = 1 ./ col_sums;
        
        % Update matrix using vectorized operations
        A = diag(R_new) * A * diag(S_new);
        
        % Update cumulative scaling vectors
        R = R .* R_new;
        S = S .* S_new';

        % Check convergence
        if all(abs(R_new - 1) <= tol) && all(abs(S_new - 1) <= tol)
            break;
        end
    end
end