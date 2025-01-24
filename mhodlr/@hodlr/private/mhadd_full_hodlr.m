function C = mhadd_full_hodlr(A, B, operator)
% Output hodlr format when inputs A, B are HODLR matrix.
% 
    [m, n] = size_t(A);
    [mb, nb] = size_t(B);
    C = A;

    % vareps = max(A.vareps, B.vareps);
    
    if isempty(A.D) & isempty(B.D) 
        if m ~= mb | n ~= nb
            error('Please enter the inputs with consistent dimensions.');
        end  

        C.A11 = mhadd_full_hodlr(A.A11, B.A11, operator);
        C.A22 = mhadd_full_hodlr(A.A22, B.A22, operator);
        
        if operator == '+'
            C.U2 = [A.U2, B.U2];
            C.V1 = [A.V1; B.V1];
            C.U1 = [A.U1, B.U1];
            C.V2 = [A.V2; B.V2];
            %[C.U2, C.V1] = A.compress(A.U2 * A.V1 + B.U2 * B.V1);
            %[C.U1, C.V2] = A.compress(A.U1 * A.V2 + B.U1 * B.V2);
        else
            C.U2 = [A.U2, -B.U2];
            C.V1 = [A.V1; B.V1];
            C.U1 = [A.U1, -B.U1];
            C.V2 = [A.V2; B.V2];
            %[C.U2, C.V1] = A.compress(A.U2 * A.V1 - B.U2 * B.V1);
            %[C.U1, C.V2] = A.compress(A.U1 * A.V2 - B.U1 * B.V2);
        end

    else
        % if operator == '+'
        %     C.D = A.D + B.D;
        %
        % else
        %    disp(A)
        %    disp(B)
        %    disp(A.D)
        %    disp(B.D)
        %    disp("-----")
        %    C.D = A.D - B.D;
        % end
        % >> % QR 
        % >> rng(0); %fix randomness
        % >> A = rand(60, 50);
        % >> depth = 88;
        % >> min_block_size = 2;
        % >> epsilon = 1e-14;
        % >> hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default
        % >> 
        % >> [Q, R] = hqr(hA, 'kressner');
        % Error using -
        % Arrays have incompatible sizes for this operation.
        if operator == '+'
            C.D = A.dense + B.dense;
        else
            C.D = A.dense - B.dense;
        end
    end 
end
