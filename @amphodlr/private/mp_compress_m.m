%% FOR TEST ONLY
function [U, V] = mp_compress_m(A, method, epsilon)
    
    if min(size(A)) == 0
        U = zeros(size(A,1),0);
        V = zeros(0,size(A,2));
        return;
    end

    global opt;
    if strcmp(method, 'svd')
        [U, S, V] = svd(full(A), 'econ');
        S = mchop(S);
        rnk = sum(abs(diag(S)) > S(1,1) * epsilon);
        U = mchop(U(:,1:rnk));
        V = S(1:rnk,1:rnk) * mchop(V(:,1:rnk)');
        
    elseif strcmp(method, 'qr')
        [U, V, P] = qr(A);
        V = mchop(V);
        rnk = sum(abs(diag(V)) > V(1,1) * epsilon);
        U = mchop(U(:, 1:rnk));
        V = mchop(V(1:rnk,:)*P');
    end
    
    U = sparse(U);
    V = sparse(V);
end

