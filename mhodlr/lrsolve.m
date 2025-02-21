function x = lrsolve(A, b, option)
  if strcmp(option, "qr")
    x = qr_solve(A, b);
  elseif strcmp(option, "lu")
    x = lu_solve(A, b);
  else
    x = pm_solve(A, b);
  end
end
