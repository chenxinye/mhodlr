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

    BL = zeros(0, 0);
    BR = zeros(0, n);
    C = zeros(0, n);
    nrm_A = hnorm(hA, 2);
    
    [Y, YBL, YBR, YC, T, A] = iter_qr(hA, BL, BR, C, nrm_A);

    % Y is HODLR matrix 
    if nargout <= 2
        Q = hdot(hdot(Y, T), Y.transpose());
        I = hodlr('eye', m, A.bottom_level, A.min_block_size);
        
        Q = sub(I, Q);
        Y = Q;
        T = A;
    end
end


function [YA, BL, YBR, YC, T, hA] = iter_qr(hA, BL, BR, C, nrm_A)
    % """hA is HODLR matrix, BL, BR, C, are dense matrices, nrm_A is a scalar."""

    % """Return: T: hodlr"""
    [m, n] = size_t(hA, 1);
    
    q = size(C, 1);
    
    if size(BL, 1) > 0
        [BL, R] = qr(BL, 0);
        BR = R*BR;
    end
    
    % disp(size(BR));
    p = size(BR, 1);
    % disp([p, m])
    if ~isempty(hA.D)
        [Y, T, R] = wyqr([hA.D; BR; C]);

        YA  = hodlr(Y(1:m, :), 0, hA.min_block_size);

        YBR = Y(m+1:m+p, :);
        YC  = Y(m+p+1:end, :);
        
        hA = hodlr(R(1:m,:), 0, hA.min_block_size);
        T = hodlr(T, 0, hA.min_block_size); % generate matrix of 0 depths
    else
        % Compute QR decomposition of first block column
        [m1, n1] = size_t(hA.A11, 1);
        BC = [BR; C];
        % disp(size(BC))
        [YA11, YBL1, YBR1, YC1, T1, hA.A11] = iter_qr(hA.A11, hA.U2, hA.V1, BC(:, 1:n1),nrm_A*hA.vareps);
        
        SL = [hdot(YA11.transpose(), hA.U1, 'dense'), YBR1', YC1'];
        SR = [hA.V2', hdot(hA.A22.transpose(), YBL1, 'dense'), BC(:,n1+1:end)'];
        
        [SL, SR] = hrank_truncate(SL, SR', nrm_A*hA.vareps);
        SR = SR';
        % disp(size(SR))

        SL = hdot(T1.transpose(), SL, 'dense');
        
        % Update second block column
        [hA.U1, hA.V2] = hrank_truncate( [hA.U1, -hdot(YA11, SL, 'dense')], [hA.V2', SR]', nrm_A*hA.vareps);
        
        hA.A22 = fusedma(hA.A22, YBL1, (-SR* ( YBR1 * SL )')', nrm_A*hA.vareps);
        BC(:, n1+1:end) = BC(:, n1+1:end) - (YC1*SL)*SR';
        
        % Compute QR decomposition of second block column
        [m2, n2] = size_t(hA.A22, 1);
        
        [YA22, YBL2, YBR2, YC2, T2, hA.A22] = iter_qr(hA.A22, zeros(0,0), zeros(0, n2), BC(:, n1+1:end), nrm_A*hA.vareps);
        
        % Set Y
        YA = hodlr;
        YA.min_block_size = hA.min_block_size;

        YA.max_level = YA.max_level;
        YA.bottom_level = max(YA11.bottom_level, YA22.bottom_level) + 1;

        YA.A11 = YA11; 
        YA.A22 = YA22; 
        YA.U2 = YBL1; 
        YA.V1 = YBR1;
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

        T.max_level = hA.max_level;
        T.bottom_level = max(T1.bottom_level, T2.bottom_level) + 1;

        T.U2 = zeros(n2, 0); 
        T.V1 = zeros(0, n1);
        T12L = [YBR1', YC1(1:p,:)', YC1(p+1:end,:)'];
        T12R = [hdot(YA22.transpose(), YBL1, 'dense'), YBR2', YC2'];
        [T12L, T12R] = hrank_truncate(T12L, T12R', hA.vareps);
        T12R = T12R';
        
        T.U1 = -hdot(T1, T12L, 'dense');
        T.V2 = hdot(T2.transpose(), T12R, 'dense')';

    end
    
end



function [Y, T, A] = wyqr(A)
    % modified from https://github.com/numpi/hm-toolbox/blob/master/%40hodlr/private/qrWY.m
    [m,n] = size(A);
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
        
        [Y1, T1, A(:, 1:n1)] = wyqr( A(:, 1:n1) );
        x = Y1*T1';
        A(:, j1:n) = A(:, j1:n) - (x*(Y1'*A(1:m, j1:n)));
        
        [Y2, T2, A(j1:m,j1:n)] = wyqr( A(j1:m,j1:n) );
        Y2 = [zeros(size(Y1,1)-size(Y2,1),size(Y2,2));Y2];
        
        Y = [Y1,Y2];
        m2 = size(T1,2);
        n2= size(T2,1);
        T3 = -x'*(Y2*T2);
        
        T = [T1 T3; zeros(n2,m2) T2];
    end
    
    end
    
function [v,beta] = house(x) 
    % from https://github.com/numpi/hm-toolbox/blob/master/%40hodlr/private/qrWY.m
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


