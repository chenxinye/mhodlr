function H = htruncate(H, rnk)
%{
    Function to truncate the off-diagonal blocks of HODLR matrix to rank ``rnk``.
    
    Parameters
    --------------------
    
    H - hodlr, mphodlr, and amphodlr
        The input of HODLR matrix.
    
    rnk - int
        the maximum rank for off-diagonal blocks.
%}
    if ~isempty(H.D)
        return
    else
        [H.U1, H.V2] = hrank_truncate_p(H.U1, H.V2, rnk);
        [H.U2, H.V1] = hrank_truncate_p(H.U2, H.V1, rnk);
    end
end 


function [Un, Vn] = hrank_truncate_p(U, V, rnk)
    [QU, RU] = qr(U, 0);
    [QV, RV] = qr(V', 0);
    
    [U,S,V] = svd(full(RU * RV'));

    [m_u, n_u] = size(U);
    [m_v, n_v] = size(V);

    min_rnk = min([m_u, n_u, m_v, n_v]);
    if min_rnk >= rnk
        min_rnk = rnk;
    end

    Un = QU * U(:,1:min_rnk) * S(1:min_rnk,1:min_rnk);
    Vn = (QV * V(:,1:min_rnk))';
end
