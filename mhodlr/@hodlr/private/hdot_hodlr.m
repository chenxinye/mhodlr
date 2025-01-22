function C = hdot_hodlr(A, B)
    if ismember(class(A), {'hodlr', 'amphodlr', 'mphodlr'}) & ismember(class(B), {'hodlr', 'amphodlr', 'mphodlr'})

        % if A.bottom_level ~= B.bottom_level
        %     error("Please ensure the two input matrices are of the same cluster structure.")
        % else 
        %     vareps = max(A.vareps, B.vareps);
        % end

        % error info: 
        %>> rng(0); %fix randomness        
        %>> A = rand(60, 50);
        %>> depth = 4;
        %>> min_block_size = 2;
        %>> epsilon = 1e-14;
        %>> hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default
        %>> 
        %>> [Q, R] = hqr(hA, 'kressner');
        %>> hdot(Q, R) % -- cause errors 
        % disp(A);
        % disp(B);
        % disp("-------");

        vareps = max(A.vareps, B.vareps);
        
        if ~isempty(A.D)
            D = A.D * B.D;
            C = hodlr(D, 0, A.min_block_size);
            C.vareps = vareps;
            [m, n] = size(D);
            C.shape = [m, n];
        else
            C = A;
            C.vareps = vareps;
            C.A11 = hadd_partial_hodlr(hdot_hodlr(A.A11, B.A11), A.U1*(A.V2*B.U2)*B.V1, '+', true);
            C.A22 = hadd_partial_hodlr(hdot_hodlr(A.A22, B.A22), A.U2*(A.V1*B.U1)*B.V2, '+', true);

            
            A12 = hdot_dense(A.A11, B.U1*B.V2) + hdot_dense(A.U1*A.V2, B.A22);
            A21 = hdot_dense(A.U2*A.V1, B.A11) + hdot_dense(A.A22, B.U2*B.V1);
            [C.U1, C.V2] = compress_m(A12, 'svd', vareps);
            [C.U2, C.V1] = compress_m(A21, 'svd', vareps);
        end

    elseif ismember(class(A), {'hodlr', 'amphodlr', 'mphodlr'})
        vareps = A.vareps;
        [mB, nB] = size(B);
        midB = ceil(nB / 2);
        
        if ~isempty(A.D)
            D = A.D * B;
            C = hodlr(D, 0, A.min_block_size);
            C.vareps = vareps;
            [m, n] = size(D);
            C.shape = [m, n];
        else
            [mA, ~, ~, ~, sv1, ~] = size_t(A, 1);
            
            C = A;
            C.vareps = vareps;
            C.A11 = hadd_partial_hodlr(hdot_hodlr(A.A11, B(1:sv1, 1:midB)), (A.U1*A.V2)*B(sv1+1:end, 1:midB), '+', true);
            C.A22 = hadd_partial_hodlr(hdot_hodlr(A.A22, B(sv1+1:end, midB+1:end)), (A.U2*A.V1)*B(1:sv1, midB+1:end), '+', true);

            A12 = hdot_dense(A.A11, B(1:sv1, midB+1:end)) + hdot_dense(A.U1*A.V2, B(sv1+1:end, midB+1:end));
            A21 = hdot_dense(A.U2*A.V1, B(1:sv1, 1:midB)) + hdot_dense(A.A22, B(sv1+1:end, 1:midB));
            [C.U1, C.V2] = compress_m(A12, 'svd', vareps);
            [C.U2, C.V1] = compress_m(A21, 'svd', vareps);
            C.shape = [mA, nB];
        end

    elseif ismember(class(B), {'hodlr', 'amphodlr', 'mphodlr'})
        vareps = B.vareps;
        
        [mA, ~] = size(A);
        midA = ceil(mA / 2);

        if ~isempty(B.D)
            D = A * B.D;
            C = hodlr(D, 0, B.min_block_size);
            C.vareps = vareps;
            [m, n] = size(D);
            C.shape = [m, n];
        else
            [~, nB, su1, ~, ~, ~] = size_t(B, 1);
            C = B;
            C.vareps = vareps;
            C.A11 = hadd_partial_hodlr(hdot_hodlr(A(1:midA, 1:su1), B.A11), A(1:midA, su1+1:end)*(B.U2*B.V1), '+', true);
            C.A22 = hadd_partial_hodlr(hdot_hodlr(A(midA+1:end, su1+1:end), B.A22), A(midA+1:end, 1:su1)*(B.U1*B.V2), '+', true);

            A12 = hdot_dense(A(1:midA, 1:su1), B.U1*B.V2) + hdot_dense(A(1:midA, su1+1:end), B.A22);
            A21 = hdot_dense(A(midA+1:end, 1:su1), B.A11) + hdot_dense(A(midA+1:end, su1+1:end), B.U2*B.V1);
            [C.U1, C.V2] = compress_m(A12, 'svd', vareps);
            [C.U2, C.V1] = compress_m(A21, 'svd', vareps);
            C.shape = [mA, nB];
        end
    else
        C = A * B;
    end
end