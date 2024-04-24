function [m, n, varargout] = hsize(H)

    if isempty(H.D) 
        su1 = size(H.U1, 1);
        su2 = size(H.U2, 1);
        sv1 = size(H.V1, 2);
        sv2 = size(H.V2, 2);
        m = su1 + su2;
        n = sv1 + sv2;
        if nargout > 2
            varargout = {su1, su2, sv1, sv2};
        end
    else
        [m, n] = size(H.D);
    end

end