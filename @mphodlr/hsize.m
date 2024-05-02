function [m, n, varargout] = hsize(H)
%% The function is to return the shape of HOLDR matrix.
%
% Parameters
% --------------------
% H - hodlr
%     Matrix in HODLR format - hodlr class.
% 
%
% Returns
% --------------------
% [m, n, su1, su2, sv1, sv2] - int
%     Indicates the size of rows and columns for H, size(H.U1, 1), size(H.U2, 1), size(H.V1, 2), size(H.V2, 2), respectively.
%     
 
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