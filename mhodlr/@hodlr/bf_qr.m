function [Q, R] = bf_qr(hA)
    if isempty(hA.D)
       [m, n, m1, m2, n1, n2]  = size_t(hA);
       
       if m ~= n
           error('Only square matrix is supported.')
       end

       H12 = hA.U1 * hA.V2;
       H21 = hA.U2 * hA.V1;
       X = hdot(H21, hA.A11.inverse)
       XTX =  hdot(X.transpose(), X)
       XXT =  hdot(X, X.transpose())
       
       IXX1 = hadd(XTX, eye(n1), '+');
       IXX2 = hadd(XXT, eye(m2), '+');

       L1 = hchol(IXX1);
       L2 = hchol(IXX2);

       R11 = hdot(L1, hA.A11);

       B1 = hadd(M12, hdot(X.transpose, hA.A22));
       R12 = lu_solve(L1, B1);

       R22 = hadd(hA.A22, hdot(X, M12), '-');
       R22 = lu_solve(L2.transpose, R22);
       
       [Q1, R11] = bf_qr(R11);
       R12 = hdot(Q1.transpose(), R12, 'dense');

       [Q2, R22] = bf_qr(R22);
       
       Q11 = hdot(L1.inverse, Q1);
       Q22 = hdot(L2.inverse, Q2);
       
       Q21 = hdot(L1.inverse, Q1, 'dense');
       Q12 = hdot(L2.inverse, Q2, 'dense');
       
       Q12 = hdot(X.transpose(), -Q12, 'dense');
       Q21 = hdot(X, Q21, 'dense');

       Q.A11 = Q11;
       Q.A22 = Q22;
       
       [Q.U1, Q.V2] = hA.compress(Q12);
       [Q.U2, Q.V1] = hA.compress(Q21);
       
       Q = hA.bottom_level;

       R.A11 = R11;
       R.A22 = R22;

       [R.U1, R.V2] = hA.compress(R12);
       
       R.U2 = zeros(m2, 1);
       R.V1 = zeros(1, m1);
    else
        [Q, R] = qr(hA.D);
        Q = hodlr(Q, 0);
        R = hodlr(R, 0);
    end
end
