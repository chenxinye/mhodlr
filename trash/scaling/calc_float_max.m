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
