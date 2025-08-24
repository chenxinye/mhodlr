function [R, S] = two_sided_diagonal_scaling(A, tol)
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