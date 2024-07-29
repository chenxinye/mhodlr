function y = kernel4(x, y)
    y = exp(-pdist2(x, y)/800);
end