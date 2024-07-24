err

x=[1e16;1e0];y=[1e16;1e7];
[(x-y)'*(x-y), sum(x.^2)-2*x'*y + sum(y.^2)]

