%% FOR TEST ONLY
function [U, V] = mp_compress_m(A, method, vareps, varargin)
    
    if nargin == 4
        norm_type = varargin{1};
    else 
        norm_type = '2';
    end 

    if min(size(A)) == 0
        U = zeros(size(A,1),0);
        V = zeros(0,size(A,2));
        return;
    end

    global opt;
    if strcmp(method, 'svd')
        [U, S, V] = svd(full(A), 'econ');
        if nargin <= 3 | norm_type == '2'
            rnk = sum(abs(diag(S)) > S(1,1) * vareps);
            S = mchop(S);
            U = mchop(U(:,1:rnk));
            V = S(1:rnk,1:rnk) * mchop(V(:,1:rnk)');
        elseif norm_type == 'fro'
            sq_dS = diag(S).^2;
            normf = sum(sq_dS);
            cusm = cumsum(sq_dS, "reverse") / normf;
            in_eq = cusm > vareps^2;
            rnk = min(sum(in_eq) + 1, size(U, 2));
            S = mchop(S);
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
    
    U = sparse(U);
    V = sparse(V);
end

