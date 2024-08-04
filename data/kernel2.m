function y = kernel2(x, y) 
    y = pdist2(x, y);
    y(find(y==0)) = 1;
    y = log(y);
end
