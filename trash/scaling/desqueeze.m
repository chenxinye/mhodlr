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