%% FOR TEST ONLY
function [U, V] = compress_m(A, method, vareps, varargin)

    if nargin == 4
        issparse = varargin{1};
        norm_type = '2';
    elseif nargin == 5
        issparse = varargin{1};
        norm_type = varargin{2};
    else 
        issparse = true;
        norm_type = '2';
    end 

    if min(size(A)) == 0
        U = zeros(size(A,1), 0);
        V = zeros(0, size(A,2));
        return;
    end

    if strcmp(method, 'svd')
        [U, S, V] = svd(full(A), 'econ');

        if nargin <= 3 | norm_type == '2'
            rnk = sum(abs(diag(S)) > S(1,1) * vareps);
            U = U(:,1:rnk);
            V = S(1:rnk,1:rnk) * V(:,1:rnk)';
        elseif norm_type == 'fro'
            sq_dS = diag(S).^2;
            normf = sum(sq_dS);
            cusm = cumsum(sq_dS, "reverse") / normf;
            in_eq = cusm > vareps^2;
            rnk = min(sum(in_eq) + 1, size(U, 2));
            U = U(:,1:rnk);
            V = S(1:rnk,1:rnk) * V(:,1:rnk)';
        else
            error("This norm type is not sopported for truncation. Please use `fro` or `2` for norm type specification.")
        end 

    elseif strcmp(method, 'rsvd')
        [U, S, V] = rsvd(A, rnk, 1)
        if nargin <= 3 | norm_type == '2'
            rnk = sum(abs(diag(S)) > S(1,1) * vareps);
            U = U(:,1:rnk);
            V = S(1:rnk,1:rnk) * V(:,1:rnk)';
        elseif norm_type == 'fro'
            sq_dS = diag(S).^2;
            normf = sum(sq_dS);
            cusm = cumsum(sq_dS, "reverse") / normf;
            in_eq = cusm > vareps^2;
            rnk = min(sum(in_eq) + 1, size(U, 2));
            U = U(:,1:rnk);
            V = S(1:rnk,1:rnk) * V(:,1:rnk)';
        else
            error("This norm type is not sopported for truncation. Please use `fro` or `2` for norm type specification.")
        end 

        
    elseif strcmp(method, 'qr')
        [U, V, P] = qr(A);
        rnk = sum(abs(diag(V)) > V(1,1) * vareps);
        U = U(:, 1:rnk);
        V = V(1:rnk,:)*P';
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
return
