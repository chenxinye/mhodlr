function C = inverse(H, varargin)
%{
    The function is to return the inverse of HOLDR matrix.
    
    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
    
    oformat - str, default = 'hodlr'
        The format of returns.
    
    algorithm - int, default=1
        The algorithm to implement inverse
    
    Returns
    --------------------
    C - hodlr | double
        Return matrix in hodlr class or dense array.
    
%}  
    
    switch nargin
        case 1 
            oformat = 'hodlr';
            algorithm = 1;

        case 2
            oformat = varargin{1};
            algorithm = 1;
            
        case 3 
            oformat = varargin{1};
            algorithm = varargin{2};
            
    end

    if algorithm == 1
        if strcmp(oformat, 'hodlr') | strcmp(oformat, 'mphodlr') | strcmp(oformat, 'amphodlr')
            C = inverse_hodlr(H);
        else
            C = inverse_dense(H);
        end
    else
        if strcmp(oformat, 'hodlr') | strcmp(oformat, 'mphodlr') | strcmp(oformat, 'amphodlr')
            C = inverse_nonrecursive_hodlr(H);
        else
            C = inverse_nonrecursive_dense(H);
        end
    end
end 

function C = inverse_nonrecursive_hodlr(H)
    %% This member method implement inverse of hodlr matrix in nonrecursive manner
    [m, n] = size_t(H);

    if m ~= n
        error('Inverse is only applied to a square HODLR matrix.');
    end

    if isempty(H.D)
        C1 = inverse_dense(H.A11);
        C2 = inverse_dense(H.A22);
        
        A12 = H.U1 * H.V2;
        A21 = H.U2 * H.V1;
        X = blkoffdiag(C1 * A12, C2* A21);
        [U, D, V] = svd(X);

        L = eye(m) - U * inv(inv(D) + V' * U) * V';
        R = blkdiag(C1, C2);
        C = L * R;
    else
        C = inv(H.D);
    end

    [md, td, mbs, ml, tp] = load_params(H, 0);
    C = hodlr(C, md, td, mbs, ml, tp);
end

function C = inverse_nonrecursive_dense(H)
    %% This member method implement inverse of hodlr matrix in nonrecursive manner
    [m, n] = size_t(H);

    if m ~= n
        error('Inverse is only applied to a square HODLR matrix.');
    end
    
    if isempty(H.D)
        C1 = inverse_nonrecursive_dense(H.A11);
        C2 = inverse_nonrecursive_dense(H.A22);
        
        A12 = H.U1 * H.V2;
        A21 = H.U2 * H.V1;
        X = blkoffdiag(C1 * A12, C2* A21);
        [U, D, V] = svd(X);

        L = eye(m) - U * inv(inv(D) + V' * U) * V';
        R = blkdiag(C1, C2);
        C = L * R;
    else
        C = inv(H.D);
    end
    
end

function C = inverse_hodlr(H)
    if ~issquare(H)
        error('Inverse is only applied to a square HODLR matrix.');
    end
    
    C = H;
    if isempty(H.D)
        X22 = inverse_hodlr(H.A22);
        A12 = H.U1 * H.V2;
        A21 = H.U2 * H.V1;
        
        C.A11  = inverse_hodlr(hadd(H.A11, hdot_dense(hdot(A12, X22), A21), '-'));
    
        [C.U1, C.V2] = compress_m(hdot_dense(hdot_dense(C.A11, -A12), X22), H.method, H.vareps, H.max_rnk, H.trun_norm_tp, H.issparse);
        C21 = -hdot_dense(hdot_dense(X22, A21), C.A11);
        [C.U2, C.V1] = compress_m(C21, H.method, H.vareps, H.max_rnk, H.trun_norm_tp, H.issparse);
        XX = hdot_dense(C21 * A12, X22);
        C.A22 = hadd(X22, XX, '-');
    else
        C.D = inv(H.D);
    end
end

function C = inverse_dense(H)
    if ~issquare(H)
        error('Inverse is only applied to a square HODLR matrix.');
    end
    
    if strcmp(class(H), 'amphodlr') | strcmp(class(H), 'hodlr') | strcmp(class(H), 'mphodlr')
        if isempty(H.D)
            X22 = inverse_dense(H.A22);
            A12 = H.U1*H.V2;
            A21 = H.U2*H.V1;
            X11 = inverse_dense(hadd_partial_hodlr(H.A11, A12 * X22 * A21 ,'-', true));
            C12 = -X11 * A12 * X22;
            C21 = -X22 * A21 * X11;
            C22 = X22 + X22 * A21 * X11 * A12 * X22;

            C = [X11, C12; C21, C22];
        else
            C = inv(H.D);
        end
    else
        C = inv(H);
    end
end
