function [x] = qr_solve(varargin)
    %{
        Compute Hx = b using LU factorization.
    
        Parameters
        --------------------
        Setting 1:
            method - str
                ''lintner'' | ''bebendorf'' | ''kressner''

            H - hodlr
                Matrix in HODLR format - hodlr class.
    
            b - double/single
    
    
        Setting 2:
            Q - hodlr
                Matrix in HODLR format - hodlr class.
    
            R - hodlr
                Matrix in HODLR format - hodlr class.
                
            b - double/single
    
        Returns
        --------------------
        x - double
            The solution.
    %}
    method = varargin{1};
    if isa(varargin{1}, "string") | isa(varargin{1}, "char")
        
        H = varargin{2}; 
        if ~(isa(H, 'hodlr') | isa(H, 'mphodlr') | isa(H, 'amphodlr'))
            error('Please ensure the first input is of a HODLR matrix.');
        end
        
        if strcmp(method, 'lintner')
            [Q, R] = lintner_qr(H);
    
        elseif strcmp(method, 'bebendorf')
            [Q, R] = bebendorf_qr(H);
    
        elseif strcmp(method, 'kressner')
            [Q, R] = kressner_qr(H);
        end
        
        b = varargin{3}; 
        y = hdot(Q.transpose(), b, 'dense');
        x = htrsu(R, y, 2);

    else
        Q = varargin{1}; 
        R = varargin{2}; 

        if ~(isa(Q, 'hodlr') | isa(Q, 'mphodlr') | isa(Q, 'mphodlr'))
            error('Please ensure the first input is of a HODLR matrix.');
        end

        if ~(isa(R, 'hodlr') | isa(R, 'mphodlr') | isa(R, 'amphodlr'))
            error('Please ensure the second input is of a HODLR matrix.');
        end
        
        b = varargin{3};
        y = hdot(Q.transpose(), b, 'dense');
        x = htrsu(U, y, 2);
    end
end