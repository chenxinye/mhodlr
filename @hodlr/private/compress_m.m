%% FOR TEST ONLY
function [U, V] = compress_m(A, method, threshold)
    if min(size(A)) == 0
        U = zeros(size(A,1),0);
        V = zeros(0,size(A,2));
        return;
    end

    if strcmp(method, 'svd')
        [U, S, V] = svd(full(A), 'econ');
        k = sum(abs(diag(S)) > S(1,1) * threshold);
        U = U(:,1:k) * S(1:k,1:k);
        V = V(:,1:k)';
    elseif strcmp(method, 'qr')
        [U, V, P] = qr(A);
        k = sum(abs(diag(V)) > V(1,1) * threshold);
        U = U(:, 1:k);
        V = V(1:k,:)*P;
    end
end