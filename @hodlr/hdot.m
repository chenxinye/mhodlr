function C = hdot(A, B, varargin) % one if left A, the other is right A
    if isa(A, 'hodlr') & isa(B, 'hodlr')
        error('Please enter correct types of inputs.');
    end

    if nargin == 2 
        oformat = 'hodlr';
    else
        if strcmp(varargin{1}, 'hodlr')
            oformat = 'hodlr';
        else
            oformat = 'matlab_default';
        end
    end

    if isa(A, 'hodlr')
        if isempty(A.D)
            [mA, nA, su1, su2, sv1, sv2] = hsize(A);
            
            [mB, nB] = size(B);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.');
            end
        
            C = zeros(mA, nB);
            y1 = hdot(A.A11, B(1:su1, :), 'double');
            y2 = A.U1 * A.V2 * B(su1+1:end, :);
            y3 = A.U2 * A.V1 * B(1:su1, :);
            y4 = hdot(A.A22, B(su1+1:end, :), 'double');
            
        
            C(1:su1, :) = y1 + y2;
            C(su1+1:end, :) = y3 + y4;
        else
            C = A.D * B;
        end

    elseif isa(B, 'hodlr')
        if isempty(B.D)
            [mB, nB, su1, su2, sv1, sv2] = hsize(B);
            
            [mA, nA] = size(A);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.')
            end
        
            C = zeros(mA, nB);
            y1 = hdot(A(:, 1:su1), B.A11, 'double');
            y2 = A(:, su1+1:end) * B.U2 * B.V1;
            y3 = A(:, 1:su1) * B.U1 * B.V2;
            y4 = hdot(A(:, su1+1:end), B.A22, 'double');
            
            C(:, 1:su1) = y1 + y2;
            C(:, su1+1:end) = y3 + y4;
        else
            C = A * B.D;
        end
        
    else
        C = A * B;
    end
    
    if strcmp(oformat, 'hodlr')
        C = hodlr(C);
    end
    
end