function A = recover(H, varargin)
%{
    The function is to recover a HODLR format into array format

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
 
    issparse - boolean
        `1` indicates returning sparse format, `0` indicates returning full arrary. 


    Returns
    --------------------
    A - double 
        Array in sparse or not.
%}
    A = recover_mat(H);
    if nargin > 1
        A = sparse(A);
    end
end

function A = recover_mat(H)
    if isempty(H.D) 
        su1 = size(H.U1, 1);
        su2 = size(H.U2, 1);
        sv1 = size(H.V1, 2);
        sv2 = size(H.V2, 2);
        rowSize = su1 + su2;
        colSize = sv1 + sv2;
        A = zeros(rowSize, colSize);
        A(1:su1, sv1+1:end) = H.U1 * H.V2;
        A(su1+1:end, 1:sv1) = H.U2 * H.V1;
        A(1:su1, 1:sv1) = recover_mat(H.A11);
        A(su1+1:end, sv1+1:end) = recover_mat(H.A22);
    else
        A = H.D;
    end
end
