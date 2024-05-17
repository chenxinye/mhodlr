%% FOR TEST ONLY
function [U, V] = compress_m(A, method, epsilon)
    if min(size(A)) == 0
        U = zeros(size(A,1),0);
        V = zeros(0,size(A,2));
        return;
    end

    if strcmp(method, 'svd')
        [U, S, V] = svd(full(A), 'econ');
        rnk = sum(abs(diag(S)) > S(1,1) * epsilon);
        U = U(:,1:rnk) * S(1:rnk,1:rnk);
        V = V(:,1:rnk)';
        
    elseif strcmp(method, 'qr')
        [U, V, P] = qr(A);
        rnk = sum(abs(diag(V)) > V(1,1) * epsilon);
        U = U(:, 1:rnk);
        V = V(1:rnk,:)*P';
    end
    
    U = sparse(U);
    V = sparse(V);
end