function y = kernel3(x, y)
    y = exp(-pdist2(x, y)/2);
end