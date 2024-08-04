function arr = prec_chain(varargin)
%{
    Build a chain for precision used, return a cell array where 
    each element contains a precision
%}
    arr = {struct};
    for i = 1:nargin
        arr{i} = varargin{i};
    end 
end