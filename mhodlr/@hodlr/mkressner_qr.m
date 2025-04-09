function [Y, T, A] = kressner_qr(hA)
    % It is important to note that all operations, except for matrix-vector
    % multiplication, are combined with low-rank truncation, as discussed above, to limit rank growth
    % in the off-diagonal blocks. 
    % [Q,R] = kressner_qr(A) computes a QR decomposition A = Q*R of a HODLR matrix
    %     such that Q is a numerically orthogonal HODLR matrix and R is an upper
    %     triangular matrix.
    %
    % [Y,T,R] = kressner_qr(A) returns the factor Q in terms of its compact WY
    %     representation: Q = I-Y*T*Y' with an upper triangular matrix T and a
    %     lower triangular matrix Y.
    %
    % We follow the paper  ``D. Kressner and A. Šusnjara. (2018). Fast QR decomposition of HODLR
    %     matrices. Technical report, September 2018.``
    % Modified from https://github.com/numpi/hm-toolbox/blob/master/%40hodlr/qr.m
    % Not that our implementation is not limited to square matrix, it can be also extended to rectangular matrix.
    [m,n] = size_t(hA);

    % if m~=n
    %     error("ValueError, please enter a square matrix. ")
    % end
    hA = hmchop(hA);

    BL = zeros(0, 0);
    BR = zeros(0, n);
    C = zeros(0, n);
    nrm_A = hnorm(hA, 2);
    
    global opt;
    set_prec(opt);

    [Y, YBL, YBR, YC, T, A] = miter_qr(hA, BL, BR, C, nrm_A);

    % Y is HODLR matrix 
    if nargout <= 2
        Q = mhdot(mhdot(Y, T), Y.transpose());
        I = hodlr('eye', m, A.bottom_level, A.min_block_size);
        
        Q = msub(I, Q);
        Y = Q;
        T = A;
    end
end


function [YA, BL, YBR, YC, T, hA] = miter_qr(hA, BL, BR, C, nrm_A)
    % """hA is HODLR matrix, BL, BR, C, are dense matrices, nrm_A is a scalar."""

    % """Return: T: hodlr"""
    [m, n] = size_t(hA, 1);
    
    q = size(C, 1);
    
    if size(BL, 1) > 0
        [BL, R] = qr(mchop(BL), 0);
        BR = mchop(mchop(R)*mchop(BR));
    end
    
    % disp(size(BR));
    p = size(BR, 1);
    % disp([p, m])
    if ~isempty(hA.D)
        [Y, T, R] = mwyqr([mchop(hA.D); BR; C]);

        Y = mchop(Y); T = mchop(T); R = mchop(R); 
        YA  = hodlr(Y(1:m, :), 0, hA.min_block_size);

        YBR = mchop(Y(m+1:m+p, :));
        YC  = mchop(Y(m+p+1:end, :));
        
        hA = hodlr(R(1:m,:), 0, hA.min_block_size);
        T = hodlr(T, 0, hA.min_block_size); % generate matrix of 0 depths
    else
        % Compute QR decomposition of first block column
        [m1, n1] = size_t(hA.A11, 1);
        BC = mchop([BR; C]);
        % disp(size(BC))
        [YA11, YBL1, YBR1, YC1, T1, hA.A11] = miter_qr(hA.A11, hA.U2, hA.V1, BC(:, 1:n1),nrm_A*hA.vareps);
        
        SL = [mhdot(YA11.transpose(), hA.U1, 'dense'), mchop(YBR1'), mchop(YC1')];
        SR = [hA.V2', mhdot(hA.A22.transpose(), YBL1, 'dense'), BC(:,n1+1:end)'];
        
        [SL, SR] = mhrank_truncate(SL, SR', nrm_A*hA.vareps);
        SR = SR';
        % disp(size(SR))

        SL = mhdot(T1.transpose(), SL, 'dense');
        SR = mchop(SR); SL = mchop(SL); 
        
        
        % Update second block column
        [hA.U1, hA.V2] = mhrank_truncate( [hA.U1, -hdot(YA11, SL, 'dense')], [hA.V2', SR]', nrm_A*hA.vareps);
        
        hA.A22 = mfusedma(hA.A22, YBL1, (-SR* ( YBR1 * SL )')', nrm_A*hA.vareps);
        BC(:, n1+1:end) = mchop(BC(:, n1+1:end)) - (YC1*SL)*SR';
        
        % Compute QR decomposition of second block column
        [m2, n2] = size_t(hA.A22, 1);
        
        [YA22, YBL2, YBR2, YC2, T2, hA.A22] = miter_qr(hA.A22, zeros(0,0), zeros(0, n2), BC(:, n1+1:end), nrm_A*hA.vareps);
        
        YC1 = mchop(YC1); YC2 = mchop(YC2); 
        % Set Y
        YA = hodlr;
        YA.min_block_size = hA.min_block_size;

        YA.max_level = YA.max_level;
        YA.bottom_level = max(YA11.bottom_level, YA22.bottom_level) + 1;

        YA.A11 = YA11; 
        YA.A22 = YA22; 
        YA.shape(1) = YA11.shape(1) + YA22.shape(1);
        YA.shape(2) = YA11.shape(2) + YA22.shape(2);
        YA.U2 = mchop(YBL1); 
        YA.V1 = mchop(YBR1);
        YA.U1 = zeros(m1,0); 
        YA.V2 = zeros(n2,0)';
        YBR = [YC1(1:p,:), YC2(1:p,:)];
        YC = [YC1(p+1:end,:), YC2(p+1:end,:)];
        
        hA.U2 = zeros(m2, 0); 
        hA.V1 = zeros(0, n1);
        
        T = hodlr; 
        T.min_block_size = hA.min_block_size;
        T.A11 = T1;
        T.A22 = T2;
        T.shape(1) = T1.shape(1) + T2.shape(1);
        T.shape(2) = T1.shape(2) + T2.shape(2);
        T.max_level = hA.max_level;
        T.bottom_level = max(T1.bottom_level, T2.bottom_level) + 1;

        T.U2 = zeros(n2, 0); 
        T.V1 = zeros(0, n1);
        T12L = [YBR1', YC1(1:p,:)', YC1(p+1:end,:)'];
        T12R = [mhdot(YA22.transpose(), YBL1, 'dense'), YBR2', YC2'];
        [T12L, T12R] = mhrank_truncate(T12L, T12R', hA.vareps);
        T12R = mchop(T12R');
        
        T.U1 = -mhdot(T1, T12L, 'dense');
        T.V2 = mhdot(T2.transpose(), T12R, 'dense')';

    end
    
end



function [Y, T, A] = mwyqr(A)
    % modified from https://github.com/numpi/hm-toolbox/blob/master/%40hodlr/private/qrWY.m
    [m,n] = size(A);
    nb = 10;

    A = mchop(A);
    
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
            [u,beta(j)] = mhouse(A(j:end,j));
            A(j:end,j:end) = mchop(A(j:end,j:end) - ( mchop(beta(j)*u) )*( mchop(u'*A(j:end,j:end) )));
            A(j+1:end,j) = 0;
            Y(j:end,j) = u;
        end

        T = zeros(n);
        for j = 1:n,
            T(1:j-1,j) = mchop(-beta(j)*T(1:j-1,1:j-1)*(Y(:,1:j-1)'*Y(:,j)));
            T(j,j) = mchop(beta(j));
        end
    else
        % Compute QR recursively a la Elmroth and Gustavson.
        n1 = floor(n/2);
        j1 = n1+1;
        
        [Y1, T1, A(:, 1:n1)] = mwyqr( A(:, 1:n1) );
        x = mchop(Y1*T1');
        A(:, j1:n) = mchop(A(:, j1:n) - mchop(x*(Y1'*A(1:m, j1:n))));
        
        [Y2, T2, A(j1:m,j1:n)] = mwyqr( A(j1:m,j1:n) );
        Y2 = mchop([zeros(size(Y1,1)-size(Y2,1),size(Y2,2));Y2]);
        
        Y = [Y1,Y2];
        m2 = size(T1,2);
        n2= size(T2,1);
        T3 = mchop(-x'*(Y2*T2));
        
        T = [T1 T3; zeros(n2,m2) T2];
    end
    
    end
    
function [v,beta] = mhouse(x) 
    % from https://github.com/numpi/hm-toolbox/blob/master/%40hodlr/private/qrWY.m
    n = length(x); 
    nrm = mchop(norm(x));
    if nrm ~= 0
        v1 = mchop(x(1)/nrm);
        v1 = mchop(v1 + sign(v1) + ( v1==0 ));
        beta = abs(v1);
        v1 = mchop(v1*nrm);
        v = mchop([1; x(2:end) / v1]);
    else
        v = mchop([1; x(2:n)]);
        beta = 0;
    end
end







