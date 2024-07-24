function y = kernel3(x, y) % for scalar
    y = repmat(x, size(x, 2), 1)' - repmat(y, size(y, 2), 1); 
    y = y.^2;
end
