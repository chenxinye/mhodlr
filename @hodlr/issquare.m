function ch = issquare(A)
    if class(A) == 'hodlr'
        su1 = size(A.U1, 1);
        su2 = size(A.U2, 1);
        sv1 = size(A.V1, 2);
        sv2 = size(A.V2, 2);
        
        m = su1 + su2;
        n = sv1 + sv2;
        ch = m == n;
        
    elseif class(A) == 'double' or 'single'
        [m, n] = size(A);
        ch = m == n;
    else
        error('Please enter correct type of input.');
    end
end