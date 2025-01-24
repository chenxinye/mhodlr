function C = sub(A, B, varargin)
    C = hadd(A, B, '-', 'hodlr');

    is_A_hodlr = ismember(class(A), {'hodlr', 'amphodlr', 'mphodlr'});
    is_B_hodlr = ismember(class(B), {'hodlr', 'amphodlr', 'mphodlr'});

    if is_A_hodlr & is_B_hodlr
        max_rnk = min(A.max_rnk, B.max_rnk);
    elseif is_A_hodlr
        max_rnk = A.max_rnk;
    else
        max_rnk = B.max_rnk;
    end
    
    if nargin > 2
        for i = 1:(nargin-2)
            C = hadd(C, varargin{i}, '-');
        end
    end
    
    C = htruncate(C, max_rnk);
end
