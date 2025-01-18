function varargout = hsize(A, varargin)
    [m, n] = size_t(A, 1);
    
    if nargin == 2
        dim = varargin{1};
    else
        dim = 0;
    end

    if dim == 1
        varargout = m;
    elseif dim == 2
        varargout = n;
    else
        varargout{1} = m;
        varargout{2} = n;
    end

end