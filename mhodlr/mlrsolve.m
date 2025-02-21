function x = lrsolve(A, b, option)
  if strcmp(option, "qr")
    x = mqr_solve(A, b);
  elseif strcmp(option, "lu")
    x = mlu_solve(A, b);
  else
    x = mpm_solve(A, b);
  end
end
