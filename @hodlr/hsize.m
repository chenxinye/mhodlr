function [m, n, varargout] = hsize(A)
    su1 = size(A.U1, 1);
    su2 = size(A.U2, 1);
    sv1 = size(A.V1, 2);
    sv2 = size(A.V2, 2);
    m = su1 + su2;
    n = sv1 + sv2;

    nout = max(nargout, 1) - 2;
    if nargout > 1
        varargout = {su1, su2, sv1, sv2};
    end

end