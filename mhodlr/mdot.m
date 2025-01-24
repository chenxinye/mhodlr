function C = mdot(A, B, varargin)
  C = mhdot(A, B, 'hodlr');
  if nargin > 2
      for i = 1:(nargin-2)
          C = mhdot(C, varargin{i}, 'hodlr');
      end
  end
end
