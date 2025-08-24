function max_value = calc_float_max(exp_bits, sig_bits)
    % Calculate the maximum representable value for a floating-point format.
    %
    % Parameters:
    %   exp_bits (int): Number of exponent bits
    %   sig_bits (int): Number of fraction (mantissa) bits
    %
    % Returns:
    %   max_value (double): Maximum representable value
    bias = 2^(exp_bits - 1) - 1;
    max_exponent = (2^exp_bits - 2) - bias;
    max_mantissa = 2 - 2^(-sig_bits);
    max_value = max_mantissa * (2^max_exponent);
end

function [R, S] = two_sided_diagonal_scaling_fp16(A, tol)
    % Two-sided diagonal scaling for a matrix A.
    %
    % Parameters:
    %   A (matrix): Input matrix
    %   tol (double): Tolerance for zero division (default: 1e-10)
    %
    % Returns:
    %   R (vector): Row scaling vector
    %   S (vector): Column scaling vector
    if nargin < 2
        tol = 1e-10;
    end
    
    [m, n] = size(A);
    
    % Compute row sums and scale
    R = sum(abs(A), 2);
    R(R < tol) = 1;
    R = 1 ./ R;
    
    % Apply row scaling
    A_scaled = A .* R;
    
    % Compute column sums and scale
    S = sum(abs(A_scaled), 1);
    S(S < tol) = 1;
    S = 1 ./ S;
end

function [R, S] = two_sided_diagonal_scaling_sym_fp16(A, tol, max_iter)
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
        max_iter = 100;
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
        A = (A .* R_new) .* S_new';
        
        % Update cumulative scaling vectors
        R = R .* R_new;
        S = S .* S_new;
        
        % Check convergence
        if all(abs(R_new - 1) <= tol) && all(abs(S_new - 1) <= tol)
            break;
        end
    end
end

function [rounded_A, params] = squeeze_fp16(A, theta)
    % Squeeze matrix A to fp16 representation.
    %
    % Parameters:
    %   A (matrix): Input matrix
    %   theta (double): Scaling factor
    %
    % Returns:
    %   rounded_A (matrix): Scaled matrix in fp16
    %   params (struct): Scaling parameters
    if nargin < 2
        theta = 1; % Default value for theta if not provided
    end
    
    params = struct();
    exp_bits = 5;
    sig_bits = 10;
    
    [params.R, params.S] = two_sided_diagonal_scaling_fp16(A);
    xmax = calc_float_max(exp_bits, sig_bits);
    
    A_tilde = diag(params.R) * A * diag(params.S);
    alpha = max(abs(A_tilde), [], 'all');
    params.mu = theta * xmax / alpha;
    
    % Convert to single (MATLAB's equivalent to fp16 for this purpose)
    rounded_A = single(params.mu * A_tilde);
end

function [rounded_A, params] = squeeze_sym_fp16(A, tol, theta)
    % Squeeze symmetric matrix A to fp16 representation.
    %
    % Parameters:
    %   A (matrix): Input symmetric matrix
    %   tol (double): Convergence tolerance (default: 0.1)
    %   theta (double): Scaling factor
    %
    % Returns:
    %   rounded_A (matrix): Scaled matrix in fp16
    %   params (struct): Scaling parameters
    if nargin < 2
        tol = 0.1;
    end
    if nargin < 3
        theta = 1; % Default value for theta if not provided
    end
    
    params = struct();
    exp_bits = 5;
    sig_bits = 10;
    
    [params.R, params.S] = two_sided_diagonal_scaling_sym_fp16(A, tol);
    xmax = calc_float_max(exp_bits, sig_bits);
    
    A_tilde = diag(params.R) * A * diag(params.S);
    alpha = max(abs(A_tilde), [], 'all');
    params.mu = theta * xmax / alpha;
    
    % Convert to single (MATLAB's equivalent to fp16 for this purpose)
    rounded_A = single(params.mu * A_tilde);
end

function A_recon = desqueeze(rounded_A, params)
    % Desqueeze matrix back to original scale.
    %
    % Parameters:
    %   rounded_A (matrix): Scaled matrix
    %   params (struct): Scaling parameters
    %
    % Returns:
    %   A_recon (matrix): Reconstructed matrix
    A_recon = diag(1 ./ params.R) * (double(rounded_A) / params.mu) * diag(1 ./ params.S);
end

% Main execution
function main()
    % Example usage
    A = randn(3, 3); % Random 3x3 matrix
    disp('Original A:');
    disp(A);
    
    % Test squeeze_fp16
    [A_rounded, params] = squeeze_fp16(A);
    A_recon = desqueeze(A_rounded, params);
    disp('Reconstructed A (squeeze_fp16):');
    disp(A_recon);
    
    % Test squeeze_sym_fp16
    A_sym = (A + A') / 2; % Ensure symmetry
    [A_rounded, params] = squeeze_sym_fp16(A_sym);
    A_recon = desqueeze(A_rounded, params);
    disp('Reconstructed A (squeeze_sym_fp16):');
    disp(A_recon);
    
    % Test GMRES
    b = randn(3, 1);
    [x1, ~] = gmres(A, b, [], 1e-5);
    disp('GMRES solution (original):');
    disp(x1);
    
    [A_rounded, params] = squeeze_fp16(A);
    [x2, ~] = gmres(A_rounded, params.mu * diag(params.R) * b, [], 1e-5);
    x2 = params.S .* x2;
    disp('GMRES solution (squeezed):');
    disp(x2);
end

% Run the main function
main();
