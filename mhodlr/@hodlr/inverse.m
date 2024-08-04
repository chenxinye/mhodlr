function C = inverse(obj, varargin)
    %% The function is to return the inverse of HOLDR matrix.
    %
    % Parameters
    % --------------------
    % H - hodlr
    %     Matrix in HODLR format - hodlr class.
    % 
    % algorithm - int, default=1
    %     The algorithm to implement inverse
    % 
    % oformat - str, default = 'hodlr'
    %     The format of returns.
    % 
    % Returns
    % --------------------
    % C - hodlr | double
    %     Return matrix in hodlr class or dense array.
    % 
    
    switch nargin
        case 1 
            algorithm = 1;
            oformat = 'hodlr';

        case 2
            algorithm = varargin{1};
            oformat = 'hodlr';

        case 3 
            algorithm = varargin{1};
            oformat = varargin{2};
    end

    if algorithm == 1
        if algorithm == 'hodlr'
            C = inverse_hodlr(obj);
        else
            C = inverse_dense(obj);
        end
    else
        if algorithm == 'hodlr'
            C = inverse_nonrecursive_hodlr(obj);
        else
            C = inverse_nonrecursive_dense(obj);
        end
    end
end 

function C = inverse_nonrecursive_hodlr(obj)
    %% This member method implement inverse of hodlr matrix in nonrecursive manner
    [m, n] = hsize(obj);

    if m ~= n
        error('Inverse is only applied to a square HODLR matrix.');
    end

    if isempty(obj.D)
        C1 = inverse_dense(obj.A11);
        C2 = inverse_dense(obj.A22);
        
        A12 = obj.U1 * obj.V2;
        A21 = obj.U2 * obj.V1;
        X = blkoffdiag(C1 * A12, C2* A21);
        [U, D, V] = svd(X);

        L = eye(m) - U * inv(inv(D) + V' * U) * V';
        R = blkdiag(C1, C2);
        C = L * R;
    else
        C = inv(obj.D);
    end

    [md, td, mbs, ml, tp] = load_params(obj, 0);
    C = hodlr(C, md, td, mbs, ml, tp);
end

function C = inverse_nonrecursive_dense(obj)
    %% This member method implement inverse of hodlr matrix in nonrecursive manner
    [m, n] = hsize(obj);

    if m ~= n
        error('Inverse is only applied to a square HODLR matrix.');
    end
    
    if isempty(obj.D)
        C1 = inverse_nonrecursive_dense(obj.A11);
        C2 = inverse_nonrecursive_dense(obj.A22);
        
        A12 = obj.U1 * obj.V2;
        A21 = obj.U2 * obj.V1;
        X = blkoffdiag(C1 * A12, C2* A21);
        [U, D, V] = svd(X);

        L = eye(m) - U * inv(inv(D) + V' * U) * V';
        R = blkdiag(C1, C2);
        C = L * R;
    else
        C = inv(obj.D);
    end
    
end

function C = inverse_hodlr(obj)
    if ~issquare(obj)
        error('Inverse is only applied to a square HODLR matrix.');
    end
    
    C = obj;
    if isempty(obj.D)
        X22 = inverse_hodlr(obj.A22);
        A12 = obj.U1 * obj.V2;
        A21 = obj.U2 * obj.V1;
        
        C.A11  = inverse_hodlr(hadd(obj.A11, hdot_dense(hdot(A12, X22), A21), '-'));
    
        [C.U1, C.V2] = compress_m(hdot_dense(hdot_dense(C.A11, -A12), X22), obj.method, obj.vareps);
        C21 = -hdot_dense(hdot_dense(X22, A21), C.A11);
        [C.U2, C.V1] = compress_m(C21, obj.method, obj.vareps);
        XX = hdot_dense(C21 * A12, X22);
        C.A22 = hadd(X22, XX, '-');
    else
        C.D = inv(obj.D);
    end
end

function C = inverse_dense(obj)
    if ~issquare(obj)
        error('Inverse is only applied to a square HODLR matrix.');
    end
    
    if strcmp(class(obj), 'amphodlr') | strcmp(class(obj), 'hodlr')
        if isempty(obj.D)
            X22 = inverse_dense(obj.A22);
            A12 = obj.U1*obj.V2;
            A21 = obj.U2*obj.V1;
            X11 = inverse_dense(hadd_partial_hodlr(obj.A11, A12 * X22 * A21 ,'-'));
            C12 = -X11 * A12 * X22;
            C21 = -X22 * A21 * X11;
            C22 = X22 + X22 * A21 * X11 * A12 * X22;

            C = [X11, C12; C21, C22];
        else
            C = inv(obj.D);
        end
    else
        C = inv(obj);
    end
end
