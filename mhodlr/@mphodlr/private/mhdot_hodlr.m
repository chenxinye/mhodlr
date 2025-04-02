function C = mhdot_hodlr(A, B)
    % Optimized serial matrix-matrix multiplication returning a HODLR object
    % A and B can be hodlr, amphodlr, mphodlr, or numeric arrays
    
    % Compute class flags once
    isA_hodlr = is_hodlr_class(A);
    isB_hodlr = is_hodlr_class(B);
    
    % Dispatch to helper function with precomputed flags
    C = mhdot_hodlr_helper(A, B, isA_hodlr, isB_hodlr);
end

function C = mhdot_hodlr_helper(A, B, isA_hodlr, isB_hodlr)
    % Helper function for optimized serial HODLR multiplication
    
    % Case 1: Both A and B are HODLR-like
    if isA_hodlr && isB_hodlr
        A = hmchop(A);
        B = hmchop(B);
        
        mA = A.shape(1);
        nB = B.shape(2);
        vareps = min(A.vareps, B.vareps);
        
        if ~isempty(A.D)
            D = mchop(A.D * B.D);
            C = mphodlr(D, 0, A.min_block_size);
            C.vareps = vareps;
            [m, n] = size(D);
            C.shape = [m, n];
        else
            % Precompute shared low-rank terms
            V2B_U2 = mchop(A.V2 * B.U2);
            V1B_U1 = mchop(A.V1 * B.U1);
            
            % Initialize C with A's structure
            C = A;
            C.vareps = vareps;
            C.shape = [mA, nB];
            
            % Recursive diagonal blocks with additions
            C.A11 = mhadd_partial_hodlr(mmhdot_hodlr_helper(A.A11, B.A11, true, true), ...
                                       A.U1 * (V2B_U2 * B.V1), '+', true);
            C.A22 = mhadd_partial_hodlr(mmhdot_hodlr_helper(A.A22, B.A22, true, true), ...
                                       A.U2 * (V1B_U1 * B.V2), '+', true);
            
            % Off-diagonal blocks with nested operations
            A12 = mhdot_dense(A.A11, B.U1 * B.V2) + A.U1 * (A.V2 * B.A22);
            A21 = mchop(A.U2 * mchop(A.V1 * B.A11)) + mhdot_dense(A.A22, B.U2 * B.V1);
            [C.U1, C.V2] = mp_compress_m(A12, 'svd', vareps);
            [C.U2, C.V1] = mp_compress_m(A21, 'svd', vareps);
        end
        
    % Case 2: A is HODLR-like, B is numeric
    elseif isA_hodlr
        A = hmchop(A);
        B = mchop(B);
        
        vareps = A.vareps;
        mA = A.shape(1);
        sv1 = size(A.V1, 2);  % Column split index
        [~, nB] = size(B);
        midB = ceil(nB / 2);
        
        if ~isempty(A.D)
            D = mchop(A.D * B);
            C = mphodlr(D, 0, A.min_block_size);
            C.vareps = vareps;
            [m, n] = size(D);
            C.shape = [m, n];
        else
            % Precompute shared terms
            V2B_mid = mchop(A.V2 * B(sv1+1:end, :));
            V1B = mchop(A.V1 * B(1:sv1, :));
            
            % Initialize C
            C = A;
            C.vareps = vareps;
            C.shape = [mA, nB];
            
            % Diagonal blocks
            C.A11 = mhadd_partial_hodlr(mhdot_hodlr_helper(A.A11, B(1:sv1, 1:midB), true, false), ...
                                       A.U1 * V2B_mid(:, 1:midB), '+', true);
            C.A22 = ,hadd_partial_hodlr(mhdot_hodlr_helper(A.A22, B(sv1+1:end, midB+1:end), true, false), ...
                                       A.U2 * V1B(:, midB+1:end), '+', true);
            
            % Off-diagonal blocks
            A12 = hdot_dense(A.A11, B(1:sv1, midB+1:end)) + A.U1 * V2B_mid(:, midB+1:end);
            A21 = A.U2 * V1B(:, 1:midB) + hdot_dense(A.A22, B(sv1+1:end, 1:midB));
            [C.U1, C.V2] = mp_compress_m(A12, 'svd', vareps);
            [C.U2, C.V1] = mp_compress_m(A21, 'svd', vareps);
        end
        
    % Case 3: A is numeric, B is HODLR-like
    elseif isB_hodlr
        
        A = mchop(A);
        B - hmchop(B);
        
        vareps = B.vareps;
        [mA, ~] = size(A);
        nB = B.shape(2);
        su1 = size(B.U1, 1);  % Row split index
        midA = ceil(mA / 2);
        
        if ~isempty(B.D)
            D = mchop(A * B.D);
            C = mphodlr(D, 0, B.min_block_size);
            C.vareps = vareps;
            [m, n] = size(D);
            C.shape = [m, n];
        else
            % Precompute shared terms
            AU2 = mchop(A(:, su1+1:end) * B.U2);
            AU1 = mchop(A(:, 1:su1) * B.U1);
            
            % Initialize C
            C = B;
            C.vareps = vareps;
            C.shape = [mA, nB];
            
            % Diagonal blocks
            C.A11 = hadd_partial_hodlr(mhdot_hodlr_helper(A(1:midA, 1:su1), B.A11, false, true), ...
                                       AU2(1:midA, :) * B.V1, '+', true);
            C.A22 = hadd_partial_hodlr(mhdot_hodlr_helper(A(midA+1:end, su1+1:end), B.A22, false, true), ...
                                       AU1(midA+1:end, :) * B.V2, '+', true);
            
            % Off-diagonal blocks
            A12 = AU1(1:midA, :) * B.V2 + hdot_dense(A(1:midA, su1+1:end), B.A22);
            A21 = hdot_dense(A(midA+1:end, 1:su1), B.A11) + AU2(midA+1:end, :) * B.V1;
            [C.U1, C.V2] = mp_compress_m(A12, 'svd', vareps);
            [C.U2, C.V1] = mp_compress_m(A21, 'svd', vareps);
        end
        
    % Case 4: Both A and B are numeric
    else
        D = mchop(A * B);
        C = mphodlr(D, 0, 2);  % Default min_block_size = 2
        [m, n] = size(D);
        C.shape = [m, n];
    end
end

function is_hodlr = is_hodlr_class(obj)
    % Helper function to check if obj is a HODLR-like class
    switch class(obj)
        case {'hodlr', 'amphodlr', 'mphodlr'}
            is_hodlr = true;
        otherwise
            is_hodlr = false;
    end
end

