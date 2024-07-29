function err = lu_err_bound(n, depth, epsilon, u, norm_H, norm_L, norm_U)
    err = (2^depth - 1) * (2 * epsilon + u) * norm_H + (2^depth - 1)*(7*epsilon + constant_gamma(n, u) +2*constant_gamma(n, u)) * norm_L * norm_U;
end



function rf = constant_gamma(n, u)
    rf = n*u / (1- n*u); 
end