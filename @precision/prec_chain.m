function arr = prec_chain(varargin)
    arr = {struct};
    for i = 1:nargin
        arr{i} = varargin{i};
    end 
end