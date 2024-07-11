function y = kernel2(x, y) % for scalar
    y = log(abs(repmat(x, size(x, 2), 1)' - repmat(y, size(y, 2), 1))); 
end
