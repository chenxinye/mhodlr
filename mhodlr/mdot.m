function C = mdot(A, B, varargin)
  C = hdot(A, B, 'hodlr');
  if nargin > 2
      for i = 1:(nargin-2)
          C = mhdot(C, varargin{i}, 'hodlr');
      end
  end
end
