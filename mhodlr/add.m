function C = add(A, B, varargin)
    C = hadd(A, B, '+', 'hodlr');
    if nargin > 2
        for i = 1:(nargin-2)
            C = hadd(C, varargin{i-2}, '+');
        end
    end
end
