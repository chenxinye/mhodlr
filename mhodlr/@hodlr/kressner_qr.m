function [Y, T, A] = kressner_qr(hA)
    % It is important to note that all operations, except for matrix-vector
    % multiplication, are combined with low-rank truncation, as discussed above, to limit rank growth
    % in the off-diagonal blocks. 
    % [Q,R] = QR(A) computes a QR decomposition A = Q*R of a HODLR matrix
    %     such that Q is a numerically orthogonal HODLR matrix and R is an upper
    %     triangular matrix.
    %
    % [Y,T,R] = QR(A) returns the factor Q in terms of its compact WY
    %     representation: Q = I-Y*T*Y' with an upper triangular matrix T and a
    %     lower triangular matrix Y.
    %
    % [1] D. Kressner and A. Šusnjara. (2018). Fast QR decomposition of HODLR
    %     matrices. Technical report, September 2018.

    [m,n] = hsize(hA);

    if m~=n
        error("ValueError, please enter a square matrix. ")
    end

    BL = zeros(0,0); 
    BR = zeros(0,n);
    C = zeros(0,n);
    nrm_A = hnorm(hA, 2);
    
    [Y, YBL, YBR, YC, T, A] = qr_iter(hA, BL, BR, C, nrm_A);

    if nargout <= 2,
        [r, c] = get_partitions(Y);
        Q = hdot(hdot(Y, T), Y.transpose())
        Q = hadd(hodlr('eye', Q.level, Q.min_block_size), Q, '-', 'hodlr')
        Y = Q;
        T = A;
    end

end

function [YA, BL, YBR, YC, T, A] = qr_iter(hA, BL, BR, C, nrm_A)
    """BL, BR, C, nrm_A are dense matrices"""
    [m, n] = hsize(hA, 1);
    q = size(C, 1);
    
    if size(BL,1 ) > 0,
        [BL,R] = qr(BL, 0);
        BR = R*BR;
    end
    
    p = size(BR, 1);
    
    if not isempty(hA.D)
        [Y, T, R] = qrWY([hA.D; BR; C]);
        YA  = hodlr(Y(1:m,:), hA.max_level, hA.min_block_size, hA.method, hA.vareps, hA.max_rnk, hA.trun_norm_tp);
        YBR = Y(m+1:m+p,:);
        YC  = Y(m+p+1:end,:);
        hA = hodlr(R(1:m,:), hA.max_level, hA.min_block_size, hA.method, hA.vareps, hA.max_rnk, hA.trun_norm_tp);
        T = hodlr(T, hA.max_level, hA.min_block_size, hA.method, hA.vareps, hA.max_rnk, hA.trun_norm_tp);
    else
        % Compute QR decomposition of first block column
        [m1, n1] = hsize(hA.A11);
        BC = [BR; C];
        [YA11, YBL1, YBR1, YC1, T1, hA.A11] = qr_iter(hA.A11, hA.U2, hA.V1, BC(:,1:n1), nrm_A);
        SL = [hdot(YA11.transpose(), hA.U1, 'dense'), YBR1', YC1'];
        SR = [hA.V2, hA.A22'*YBL1, BC(:,n1+1:end)'];
        [SL, SR] = hrank_truncate(SL, SR, nrm_A);
        SL = T1'*SL;
        
        % Update second block column
        [hA.U1, hA.V2] = hrank_truncate( [hA.U1, -YA11*SL] , [hA.V2, SR], nrm_A );
        hA.A22 = fusedma(A.A22, YBL1, -SR* ( YBR1 * SL )', nrm_A);
        BC(:, n1+1:end) = BC(:, n1+1:end) - (YC1*SL)*SR';
        
        % Compute QR decomposition of second block column
        [m2, n2] = size(hA.A22);
        
        [YA22, YBL2, YBR2, YC2, T2, hA.A22] = qr_iter(hA.A22, zeros(0,0), zeros(0,n2), BC(:,n1+1:end), nrm_A);
        
        % Set Y
        YA = hodlr;
        YA.A11 = YA11; 
        YA.A22 = YA22; 
        YA.sz = hsize(hA);
        YA.U2 = YBL1; 
        YA.V1 = YBR1';
        YA.U1 = zeros(m1,0); 
        YA.V2 = zeros(n2,0);
        
        YBR = [YC1(1:p,:), YC2(1:p,:)];
        YC = [YC1(p+1:end,:), YC2(p+1:end,:)];
        
        % Clean up A
        hA.U2 = zeros(m2,0); 
        hA.V1 = zeros(n1,0);
        
        % Set T
        T = hodlr; 
        T.sz = [n1+n2, n1+n2];
        T.A11 = T1; 
        T.A22 = T2;
        T.U21 = zeros(n2,0); T.V21 = zeros(n1,0);
        T12L = [YBR1', YC1(1:p,:)', YC1(p+1:end,:)'];
        T12R = [YA22'*YBL1, YBR2', YC2'];
        [T12L, T12R] = hrank_truncate(T12L, T12R, 1);
        T.U12 = -T1*T12L;
        T.V12 = T2'*T12R;
    end
    
    end



function [Y, T, A] = qrWY(A)
    [m,n] = size(A);
    if m < n,
        error('Input matrix must have more rows than columns.');
    end
    
    nb = 10;
    
    if n <= nb,
        Y = zeros(m,n);
        beta = zeros(n,1);
        if (m == n)
            beta(end) = 0;
            Y(m,n)=1;
            n1  = n-1;
        else
            n1 = n;
        end
        for j = 1:n1,
            [u,beta(j)] = house(A(j:end,j));
            A(j:end,j:end) = A(j:end,j:end) - ( beta(j)*u )*( u'*A(j:end,j:end) );
            A(j+1:end,j) = 0;
            Y(j:end,j) = u;
        end
        T = zeros(n);
        for j = 1:n,
            T(1:j-1,j) = -beta(j)*T(1:j-1,1:j-1)*(Y(:,1:j-1)'*Y(:,j));
            T(j,j) = beta(j);
        end
    else
        % Compute QR recursively a la Elmroth and Gustavson.
        n1 = floor(n/2);
        j1 = n1+1;
        
        [Y1, T1, A(:, 1:n1)] = qrWY( A(:, 1:n1) );
        x = Y1*T1';
        A(:, j1:n) = A(:, j1:n) - (x*(Y1'*A(1:m, j1:n)));
        
        [Y2, T2, A(j1:m,j1:n)] = qrWY( A(j1:m,j1:n) );
        Y2 = [zeros(size(Y1,1)-size(Y2,1),size(Y2,2));Y2];
        
        Y = [Y1,Y2];
        m2 = size(T1,2);
        n2= size(T2,1);
        T3 = -x'*(Y2*T2);
        
        T = [T1 T3; zeros(n2,m2) T2];
    end
    
    end
    
    function [v,beta] = house(x)
    % HOUSE
    %
    % Given a vector x in R^n, this computes v in R^n
    % and beta in R, such that (eye - beta v*v') x = [* 0]';
    %
    
    n = length(x);
    nrm = norm(x);
    if nrm ~= 0
        v1 = x(1)/nrm;
        v1 = v1 + sign(v1) + ( v1==0 );
        beta = abs(v1);
        v1 = v1*nrm;
        v = [1; x(2:end) / v1];
    else
        v = [1; x(2:n)];
        beta = 0;
    end
    
    end



    function [p, q] = get_partitions(hA)
       
        if ~isempty(hA.D)
            [m, n] = hsize(hA);
        else
            [m1, n1] = get_partitions(hA.A11);
            [m2, n2] = get_partitions(hA.A22);
            
            m = max(length(m1), length(m2));
            
            m1 = [ m1, ones(1, m - length(m1)) * m1(end) ];
            m2 = [ m2, ones(1, m - length(m2)) * m2(end) ];
            n1 = [ n1, ones(1, m - length(n1)) * n1(end) ];
            n2 = [ n2, ones(1, m - length(n2)) * n2(end) ];
            
            m = [ m1, m2 + m1(end) ];
            n = [ n1, n2 + n1(end) ];
        end
        
    end