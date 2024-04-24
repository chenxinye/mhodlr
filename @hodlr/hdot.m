function C = hdot(A, B, varargin) 
    if nargin == 2 
        oformat = 'hodlr';
    else
        if strcmp(varargin{1}, 'hodlr')
            oformat = 'hodlr';
        else
            oformat = 'double';
        end
    end

    C = hdot_double(A, B); 
    
    if strcmp(oformat, 'hodlr')
        C = hodlr(C);
    end
    
end