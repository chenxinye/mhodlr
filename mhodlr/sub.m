function C = sub(A, B, varargin)
    C = hadd(A, B, '-', 'hodlr');

    max_rnk = min(A.max_rnk, B.max_rnk);
    
    if nargin > 2
        for i = 1:(nargin-2)
            C = hadd(C, varargin{i-2}, '-');
        end
    end
    
    C = htruncate(C, max_rnk);
end
