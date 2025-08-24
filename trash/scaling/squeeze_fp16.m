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
    
    [params.R, params.S] = two_sided_diagonal_scaling(A);
    xmax = calc_float_max(exp_bits, sig_bits);
    
    A_tilde = diag(params.R) * A * diag(params.S);
    alpha = max(abs(A_tilde), [], 'all');
    params.mu = theta * xmax / alpha;
    
    % Convert to single (MATLAB's equivalent to fp16 for this purpose)
    rounded_A = single(params.mu * A_tilde);
end