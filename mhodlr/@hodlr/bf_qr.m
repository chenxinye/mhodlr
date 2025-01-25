function [Q, R] = bf_qr(hA)
    if isempty(hA.D)
       [m, n, m1, m2, n1, n2]  = size_t(hA);
       
       if m ~= n
           error('Only square matrix is supported.')
       end

       M12 = hA.U1 * hA.V2;
       M21 = hA.U2 * hA.V1;
       X = hdot(M21, hA.A11.inverse())
       XTX =  hdot(X.transpose(), X)
       XXT =  hdot(X, X.transpose())
       disp(size(XTX.dense))

       disp(size(XXT.dense))
       IXX1 = hadd(XTX, eye(n1), '+');
       IXX2 = hadd(XXT, eye(m2), '+');

       L1 = hchol(IXX1);
       L2 = hchol(IXX2);

       R11 = hdot(L1, hA.A11);

       R12 = M12 + X.transpose().dense * hA.A22.dense;
       R12 = hdot(L1.transpose().inverse(), R12, 'dense');

       R22 = hadd(hA.A22, hdot(X, M12), '-');
       R22 = hdot(L2.transpose().inverse.dense, R22);
       
       [Q1, R11] = bf_qr(R11);
       R12 = hdot(Q1.transpose(), R12, 'dense');

       [Q2, R22] = bf_qr(R22);
       
       Q11 = hdot(L1.inverse(), Q1, 'dense');
       Q12 = hdot(L2.inverse(), Q2, 'dense');
       Q21 = hdot(L1.inverse(), Q1, 'dense');
       Q22 = hdot(L2.inverse(), Q2, 'dense');

       Q12 = hdot(X.transpose(), -Q12, 'dense');
       Q21 = hdot(X, Q21, 'dense');
        
       Q = hodlr([Q11, Q12; Q21, Q22], hA.bottom_level);
       R = hodlr([recover(R11), R12; zeros(m2, m1), recover(R22)], hA.bottom_level);

    else
        [Q, R] = qr(hA.D);
        Q = hodlr(Q, 0);
        R = hodlr(R, 0);
    end
end
