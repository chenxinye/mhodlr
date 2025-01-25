function [Q, R] = bebendorf_qr(hA)
    if isempty(hA.D)
       [m, n, m1, m2, n1, n2]  = size_t(hA);
       
       if m ~= n
           error('Only square matrix is supported.')
       end

       M12 = hA.U1 * hA.V2;
       M21 = hA.U2 * hA.V1;
       X = hdot(M21, hA.A11.inverse()); 
  
       IXX1 = hadd(eye(n1), hdot(X.transpose(), X), '+');
       IXX2 = hadd(eye(m2), hdot(X, X.transpose()), '+');

       L1 = hchol(IXX1);
       L2 = hchol(IXX2);

       R11 = hdot(L1, hA.A11);
       R12 = hadd(M12, hdot(X.transpose(), hA.A22), '+');
       R12 = hdot(L1.transpose().inverse(), R12);

       R22 = hadd(hA.A22, hdot(X, M12), '-');
       R22 = hdot(L2.transpose().inverse(), R22);
       
       disp(R11);
       [Q1, R11] = bebendorf_qr(R11);
       [Q2, R22] = bebendorf_qr(R22);
       
       Q11 = hdot(L1.inverse(), Q1);
       Q22 = hdot(L2.inverse(), Q2);

       Q12 = hdot(L2.inverse(), Q2, 'dense');
       Q21 = hdot(L1.inverse(), Q1, 'dense');
       
       Q12 = hdot(X.transpose(), -Q12, 'dense');
       Q21 = hdot(X, Q21, 'dense');
        
       Q = hA;
       R = hA;
       
       Q.A11 = Q11; 
       Q.A22 = Q22;
       [Q.U1, Q.V1] = hA.compress(Q12);
       [Q.U2, Q.V2] = hA.compress(Q21);

       R.A11 = R11; 
       R.A22 = R22;

       R.U2 = zeros(m2, 1);
       R.V1 = zeros(1, m1);

       R.U1 = Q1.transpose().dense;
       R.V2 = R12;
       
       Q.bottom_level = hA.bottom_level;
       R.bottom_level = hA.bottom_level;
       
    else
        [Q, R] = qr(hA.D);
        Q = hodlr(Q, 0);
        R = hodlr(R, 0);
    end
end
