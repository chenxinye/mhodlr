function A = recover(hodlrA)
    A = recover_mat(hodlrA);
end

function A = recover_mat(hodlrA)
    if isempty(hodlrA.D) % sum(size(hodlrA.D)) == 0
        su1 = size(hodlrA.U1, 1);
        su2 = size(hodlrA.U2, 1);
        sv1 = size(hodlrA.V1, 2);
        sv2 = size(hodlrA.V2, 2);
        rowSize = su1 + su2;
        colSize = sv1 + sv2;
        A = zeros(rowSize, colSize);
        A(1:su1, sv1+1:end) = hodlrA.U1 * hodlrA.V2;
        A(su1+1:end, 1:sv1) = hodlrA.U2 * hodlrA.V1;
        A(1:su1, 1:sv1) = recover_mat(hodlrA.A11);
        A(su1+1:end, sv1+1:end) = recover_mat(hodlrA.A22);
    else
        A = hodlrA.D;
    end
end