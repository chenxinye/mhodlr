function C = mhdot_dense(A, B)
    % Optimized matrix-matrix multiplication for HODLR and numeric matrices
    % A and B can be hodlr, amphodlr, mphodlr, or numeric arrays
    
    % Compute class flags once
    isA_hodlr = is_hodlr_class(A);
    isB_hodlr = is_hodlr_class(B);
    
    % Dispatch to helper function with precomputed flags
    C = mhdot_dense_helper(A, B, isA_hodlr, isB_hodlr);
end

function C = mhdot_dense_helper(A, B, isA_hodlr, isB_hodlr)
    % Helper function with precomputed HODLR flags
    
    % Case 1: Both A and B are HODLR-object
    if isA_hodlr && isB_hodlr
        A = hmchop(A);
        B = hmchop(B);
        [mA, nA] = A.shape;
        [mB, nB] = B.shape;
        % [mA, nA] = size_t(A);
        % [mB, nB] = size_t(B);
        
        if nA ~= mB
            error('Please enter the inputs with correct dimensions.');
        end
        
        if isempty(A.D) && isempty(B.D)
            % Precompute low-rank products once
            A12_V2B = mchop(A.V2 * B.U2);
            A21_V1B = mchop(A.V1 * B.U1);
            
            % Recursive calls with nested operations
            C11 = mhdot_dense_helper(A.A11, B.A11, true, true) + A.U1 * (A12_V2B * B.V1);
            C12 = mhdot_dense_helper(A.A11, mchop(B.U1 * B.V2), true, false) + mchop(A.U1 * mchop(A.V2 * B.A22));
            C21 = mhdot_dense_helper(mchop(A.U2 * A.V1), B.A11, false, true) + mchop(A.U2 * mchop(A21_V1B * B.V2));
            C22 = A.U2 * (A.V1 * B.U1 * B.V2) + mhdot_dense_helper(A.A22, B.A22, true, true);
            
            C = [C11, C12; C21, C22];

        elseif ~isempty(A.D)
            C = mhdot_dense_helper(A.D, B, false, isB_hodlr);
        else
            C = mhdot_dense_helper(A, B.D, isA_hodlr, false);
        end
        
    % Case 2: A is HODLR-like, B is numeric
    elseif isA_hodlr
        A = hmchop(A);

        if isempty(A.D)
            [mA, nA, su1, su2, sv1, sv2] = size_t(A);
            [mB, nB] = size(B);
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.');
            end
            
            % Preallocate C and compute in-place
            C = zeros(mA, nB);
            V2B = mchop(A.V2 * mchop(B(sv1+1:end, :)));
            V1B = mchop(A.V1 * mchop(B(1:sv1, :)));
            C(1:su1, :) = mhdot_dense_helper(A.A11, B(1:sv1, :), true, false) + mchop(A.U1 * mchop(V2B));
            C(su1+1:end, :) = mchop(A.U2 * V1B) + mhdot_dense_helper(A.A22, B(sv1+1:end, :), true, false);
        else
            C = mchop(A.D * B);
        end
        
    % Case 3: A is numeric, B is HODLR-like
    elseif isB_hodlr
        B = hmchop(B);

        if isempty(B.D)
            [mB, nB, su1, su2, sv1, sv2] = size_t(B);
            [mA, nA] = size(A);
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.');
            end
            
            % Preallocate C and compute in-place
            C = zeros(mA, nB);
            AU1 = mchop(mchop(A(:, 1:su1)) * B.U1);
            AU2 = mchop(mchop(A(:, su1+1:end)) * B.U2);
            C(:, 1:sv1) = mhdot_dense_helper(A(:, 1:su1), B.A11, false, true) + mchop(mchop(AU2) * B.V1);
            C(:, sv1+1:end) = mchop(mchop(AU1) * B.V2) + mhdot_dense_helper(A(:, su1+1:end), B.A22, false, true);
        else
            C = mchop(A * B.D);
        end
        
    % Case 4: Both A and B are numeric
    else
        C = mchop(mchop(A) * mchop(B));
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
