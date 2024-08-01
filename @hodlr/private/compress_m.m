%% FOR TEST ONLY
function [U, V] = compress_m(A, method, vareps, varargin)

    if nargin == 4
        norm_type = varargin{4};
    else 
        norm_type = '2';
    end 

    if min(size(A)) == 0
        U = zeros(size(A,1), 0);
        V = zeros(0, size(A,2));
        return;
    end

    if strcmp(method, 'svd')
        if nargin <= 3 | norm_type == '2'
            [U, S, V] = svd(full(A), 'econ');
            rnk = sum(abs(diag(S)) > S(1,1) * vareps);
            U = U(:,1:rnk);
            V = S(1:rnk,1:rnk) * V(:,1:rnk)';
        elseif norm_type == 'fro'
            [U, S, V] = svd(full(A), 'econ');
            normf = sum(diag(S).^2);
            cusm = (normf - cumsum(diag(S).^2)) / normf;
            rnk = sum(cusm > vareps^2) + 1;
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
    
    U = sparse(U);
    V = sparse(V);
end