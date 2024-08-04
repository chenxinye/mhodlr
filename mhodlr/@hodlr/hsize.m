function varargout = hsize(H, varargin)
%{
    The function is to return the shape of HOLDR matrix.

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
    
    oformat - int, default=1
        If input is not the leafnode of the HODLR matrix:
            1: Outputs of varying number [m, n, m1, m2, n1, n2] for hierarchical block matrices.
            2: Outputs of varying number [m1, m2, n1, n2] for hierarchical block matrices.
        
        Otherwise:
            Outputs of varying number [m, n] for hierarchical block matrices.

    Returns
    --------------------
    [m, n, su1, su2, sv1, sv2] - int
        Indicates the size of rows and columns for H, size(H.U1, 1), size(H.U2, 1), size(H.V1, 2), size(H.V2, 2), respectively.
%}     
    if ~(isa(H, 'hodlr') | isa(H, 'mphodlr'))
        error('Please ensure the first input is of a HODLR matrix.');
    end
    
    if isempty(H.D) 
        su1 = size(H.U1, 1);
        su2 = size(H.U2, 1);
        sv1 = size(H.V1, 2);
        sv2 = size(H.V2, 2);
        m = su1 + su2;
        n = sv1 + sv2;

        if nargin > 1
            otype = varargin{1};
        else
            otype = 1;
        end
        
        if otype == 1
            if nargout == 1
                varargout = {m};
    
            elseif nargout == 2
                varargout = {m, n};
    
            elseif nargout == 3
                varargout = {m, n, su1};

            elseif nargout == 4
                varargout = {m, n, su1, su2};

            elseif nargout == 5
                varargout = {m, n, su1, su2, sv1};

            else
                varargout = {m, n, su1, su2, sv1, sv2};
            end
        else
            if nargout == 1
                varargout = {su1};
    
            elseif nargout == 2
                varargout = {su1, su2};
    
            elseif nargout == 3
                varargout = {su1, su2, sv1};

            else
                varargout = {su1, su2, sv1, sv2};
            end
        end

    else
        [m, n] = size(H.D);

        if nargout == 1
            varargout = {m};

        elseif nargout == 2
            varargout = {m, n};
        end
    end

end