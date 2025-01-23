function C = minverse(H, varargin)
%{
    The function is to return the minverse of HOLDR matrix.
    
    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
    
    prec - precision
        Precision to simulate the factorization.

    oformat - str, default = 'hodlr'
        The format of returns.
        
    algorithm - int, default=1
        The algorithm to implement HODLR inverse
    
    
    Returns
    --------------------
    C - hodlr | double
        Return matrix in hodlr class or dense array.
    
%}     
    global opt;
    set_prec(opt);

    if nargin ==  1
        oformat = 'hodlr';
        algorithm = 1;
            
    elseif nargin ==  2
        oformat = varargin{1};
        algorithm = 1;
            
    elseif nargin > 2
        oformat = varargin{1};
        algorithm = varargin{2};
    end

    
    if algorithm == 1
        if strcmp(oformat, 'hodlr') | strcmp(oformat, 'mphodlr') | strcmp(oformat, 'amphodlr')
            C = minverse_hodlr(H);
        else
            C = minverse_dense(H);
        end
    else
        if strcmp(oformat, 'hodlr') | strcmp(oformat, 'mphodlr') | strcmp(oformat, 'amphodlr')
            C = minverse_nonrecursive_hodlr(H);
        else
            C = minverse_nonrecursive_dense(H);
        end
    end
end 

function C = minverse_nonrecursive_hodlr(H)
    %% This member method implement minverse of hodlr matrix in nonrecursive manner
    [m, n] = size_t(H);

    if m ~= n
        error('minverse is only applied to a square HODLR matrix.');
    end

    if isempty(H.D)
        C1 = mchop(minverse_dense(H.A11));
        C2 = mchop(minverse_dense(H.A22));
        
        A12 = mchop(mchop(H.U1) * mchop(H.V2));
        A21 = mchop(mchop(H.U2) * mchop(H.V1));
        X = mchop(blkoffdiag(C1 * A12, C2* A21));
        [U, D, V] = svd(mchop(X));

        U = mchop(U);
        D = mchop(D);
        V = mchop(V);

        L = eye(m) - U * inv(inv(D) + V' * U) * V';
        R = blkdiag(C1, C2);
    
        L = mchop(L);
        R = mchop(R);
        C = L * R;
        C = mchop(C);

    else
        C = inv(mchop(H.D));
        C = mchop(C);
    end

    [md, td, mbs, ml, tp] = load_params(H, 0);
    C = hodlr(C, md, td, mbs, ml, tp);
end

function C = minverse_nonrecursive_dense(H)
    %% This member method implement minverse of hodlr matrix in nonrecursive manner
    [m, n] = size_t(H);

    if m ~= n
        error('minverse is only applied to a square HODLR matrix.');
    end

    if isempty(H.D)
        C1 = mchop(minverse_nonrecursive_dense(H.A11));
        C2 = mchop(minverse_nonrecursive_dense(H.A22));
        
        A12 = mchop(mchop(H.U1) * mchop(H.V2));
        A21 = mchop(mchop(H.U2) * mchop(H.V1));
        X = mchop(blkoffdiag(C1 * A12, C2* A21));
        [U, D, V] = svd(mchop(X));

        U = mchop(U);
        D = mchop(D);
        V = mchop(V);

        L = eye(m) - U * inv(inv(D) + V' * U) * V';
        R = blkdiag(C1, C2);
        L = mchop(L);
        R = mchop(R);
        C = L * R;
        C = mchop(C);
    else
        C = mchop(inv(mchop(H.D)));
    end
    
end

function C = minverse_hodlr(H)
    if ~issquare(H)
        error('minverse is only applied to a square HODLR matrix.');
    end
    
    C = H;
    if isempty(H.D)
        X22 = minverse_hodlr(H.A22);
        A12 = mchop(mchop(H.U1) * mchop(H.V2));
        A21 = mchop(mchop(H.U2) * mchop(H.V1));
        
        C.A11  = minverse_hodlr(hadd(H.A11, mhdot_dense(hdot(A12, X22), A21), '-'));
    
        [C.U1, C.V2] = compress_m(mhdot_dense(mhdot_dense(C.A11, -A12), X22), H.method, H.vareps, H.max_rnk, H.trun_norm_tp, H.issparse);
        C.U1 = mchop(C.U1);
        C.V2 = mchop(C.V2);

        C21 = mchop(-mhdot_dense(mhdot_dense(X22, A21), C.A11));
        [C.U2, C.V1] = compress_m(C21, H.method, H.vareps, H.max_rnk, H.trun_norm_tp, H.issparse);
        C.U2 = mchop(C.U2);
        C.V1 = mchop(C.V1);
        XX = mchop(mhdot_dense(mchop(mchop(C21) * mchop(A12)), X22));
        C.A22 = hadd(X22, XX, '-');
    else
        C.D = mchop(inv(mchop(H.D)));
    end
end

function C = minverse_dense(H)
    if ~issquare(H)
        error('minverse is only applied to a square HODLR matrix.');
    end

    if strcmp(class(H), 'amphodlr') | strcmp(class(H), 'hodlr') | strcmp(class(H), 'mphodlr')
        if isempty(H.D)
            X22 = minverse_dense(H.A22);
            A12 = mchop(mchop(H.U1)*mchop(H.V2));
            A21 = mchop(mchop(H.U2)*mchop(H.V1));
            X11 = mchop(minverse_dense(hadd_partial_hodlr(H.A11, mchop(mchop(A12 * X22) * A21) ,'-', true)));
            C12 = mchop(mchop(-X11 * A12) * X22);
            C21 = mchop(mchop(-X22 * A21) * X11);
            C22 = mchop(X22 + mchop(mchop(mchop(mchop(X22 * A21) * X11) * A12) * X22));

            C = [X11, C12; C21, C22];
        else
            C = mchop(inv(mchop(H.D)));
        end
    else
        C = mchop(inv(mchop(H)));
    end
end
        