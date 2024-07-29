function err = lu_err_bound2(depth, epsilon, norm_H, norm_L, norm_U)
    err = 2 * (2^depth - 1) * epsilon * norm_H + 11 * (2^depth - 1) * epsilon * norm_L * norm_U;
end
