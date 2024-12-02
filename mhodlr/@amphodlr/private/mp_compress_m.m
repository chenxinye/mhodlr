%% FOR TEST ONLY
function [U, V] = mp_compress_m(A, method, vareps, varargin)
    
    if nargin == 4
        max_rnk = varargin{1};
        norm_type = '2';
        issparse = true;
    elseif nargin == 5
        max_rnk = varargin{1};
        norm_type = varargin{2};
        issparse = true;
    elseif nargin == 6
        max_rnk = varargin{1};
        norm_type = varargin{2};
        issparse = varargin{3};
    else 
        max_rnk = 999;
        norm_type = '2';
        issparse = true;
    end 

    [m, n] = size(A);
    min_size = min([m, n]);

    if max_rnk > min_size
        max_rnk = min_size
    end

    if min_size == 0
        U = zeros(m, 0);
        V = zeros(0, n);
        return;
    end

    global opt;
    if strcmp(method, 'svd')
        [U, S, V] = svd(full(A), 'econ');
        if nargin <= 3 | norm_type == '2'
            rnk = sum(abs(diag(S)) > S(1,1) * vareps);
            if rnk > max_rnk
                rnk = max_rnk;
            end
            U = mchop(U(:,1:rnk));
            V = S(1:rnk,1:rnk) * mchop(V(:,1:rnk)');
        elseif norm_type == 'fro'
            sq_dS = diag(S).^2;
            normf = sum(sq_dS);
            cusm = cumsum(sq_dS, "reverse") / normf;
            in_eq = cusm > vareps^2;
            rnk = min(sum(in_eq) + 1, size(U, 2));
            if rnk > max_rnk
                rnk = max_rnk;
            end
            U = mchop(U(:,1:rnk));
            V = S(1:rnk,1:rnk) * mchop(V(:,1:rnk)');
        else
            error("This norm type is not sopported for truncation. Please use `fro` or `2` for norm type specification.")
        end 

    elseif strcmp(method, 'rsvd')
        [U, S, V] = rsvd(A, max_rnk, 1)
        if nargin <= 3 | norm_type == '2'
            rnk = sum(abs(diag(S)) > S(1,1) * vareps);
            U = mchop(U(:,1:rnk));
            V = S(1:rnk,1:rnk) * mchop(V(:,1:rnk)');
        elseif norm_type == 'fro'
            sq_dS = diag(S).^2;
            normf = sum(sq_dS);
            cusm = cumsum(sq_dS, "reverse") / normf;
            in_eq = cusm > vareps^2;
            rnk = min(sum(in_eq) + 1, size(U, 2));
            U = mchop(U(:,1:rnk));
            V = S(1:rnk,1:rnk) * mchop(V(:,1:rnk)');
        else
            error("This norm type is not sopported for truncation. Please use `fro` or `2` for norm type specification.")
        end 

    elseif strcmp(method, 'qr')
        [U, V, P] = qr(A);
        V = mchop(V);
        rnk = sum(abs(diag(V)) > V(1,1) * vareps);
        U = mchop(U(:, 1:rnk));
        V = mchop(V(1:rnk,:)*P');
    end
    
    if issparse
        U = sparse(U);
        V = sparse(V);
    end
end





function [U, S, V] = rsvd(A, k, p)
    n = size(A, 2);
    Omega = randn(n, k+p);
    [Q,~,~] = qr(A*Omega, "econ");
    [U,S,V]  = svd(Q'*A);
    U = Q*U(:,1:k);
    S = S(1:k,1:k);
    V = V(:,1:k);
end
