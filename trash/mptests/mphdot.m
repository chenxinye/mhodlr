function C = mphdot(A, B, prec_settings, varargin) 
    %{
        Compute dot product of A and B in mixed precision.
    
        Parameters
        --------------------
        A - hodlr | double
            Matrix in HODLR format or dense array.
        
        B - hodlr | double
            Matrix in HODLR format or dense array.
    
        oformat - str, default='hodlr'
            Output format: 'hodlr' or 'dense'.
    
        Returns
        --------------------
        C - hodlr | double
            The matrix of product. 
    %}
     
        if nargin == 3
            oformat = 'hodlr';
        else
            if strcmp(varargin{1}, 'hodlr')
                oformat = 'hodlr';
            else
                oformat = 'dense';
            end
        end
        
        if strcmp(class(A), 'hodlr') | strcmp(class(A), 'mphodlr') 
            max_depth_A = A.max_level;
        else
            max_depth_A = 0;
        end

        if strcmp(class(B), 'hodlr') | strcmp(class(B), 'mphodlr') 
            max_depth_B = B.max_level;
        else
            max_depth_B = 0;
        end

        max_depth = max(max_depth_A, max_depth_B);
        
        if max_depth - 1 > length(prec_settings)
            warning(['The number of precisions used are less than ' ...
                    'the maximum tree level that can achieve. The remaining' ...
                    ' level will use the working precision for simulation.']); 
            for i = length(prec_settings):max_depth - 1
                prec_settings{i} = precision('d');
            end
        end

        C = mphdot_dense(A, B, prec_settings); 
        
        if strcmp(oformat, 'hodlr')
            C = hodlr(C);
        end
        
    end


function C = mphdot_dense(A, B, prec_settings)
    if isa(A, 'mphodlr') & isa(B, 'mphodlr')
        [mA, nA] = hsize(A);
        [mB, nB] = hsize(B);
        if nA ~= mB
            error('Please enter the inputs with correct dimensions.');
        end

        disp(A.level);
        set_prec(prec_settings{A.level});
        
        if isempty(A.D) & isempty(B.D)
            A12 = mchop(mchop(A.U1) * mchop(A.V2));
            A21 = mchop(mchop(A.U2) * mchop(A.V1));
            
            B12 = mchop(mchop(B.U1) * mchop(B.V2));
            B21 = mchop(mchop(B.U2) * mchop(B.V1));
                
            C11 = mphdot_dense(A.A11, B.A11, prec_settings) + mchop(A12 * B21);
            C12 = mphdot_dense(A.A11, B12, prec_settings) + mphdot_dense(A12, B.A22, prec_settings);
            C21 = mphdot_dense(A21, B.A11, prec_settings) + mphdot_dense(A.A22, B21, prec_settings);
            C22 = mchop(A21 * B12) + mphdot_dense(A.A22, B.A22, prec_settings);

            C = [C11, C12; C21, C22];
        
        elseif ~isempty(A.D)
            C = mphdot_dense(A.D, B, prec_settings);
        else
            C = mphdot_dense(A, B.D, prec_settings);
        end
            
    elseif isa(A, 'mphodlr')
        set_prec(prec_settings{A.level});

        if isempty(A.D)
            [mA, nA, su1, su2, sv1, sv2] = hsize(A);
            [mB, nB] = size(B);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.');
            end
        
            C = zeros(mA, nB);
            y1 = mphdot_dense(A.A11, B(1:sv1, :), prec_settings);
            y2 = mchop(mchop(mchop(A.U1) * mchop(A.V2)) * mchop(B(sv1+1:end, :)));
            y3 = mchop(mchop(mchop(A.U2) * mchop(A.V1)) * mchop(B(1:sv1, :)));
            y4 = mphdot_dense(A.A22, B(sv1+1:end, :), prec_settings);
            
            C(1:su1, :) = mchop(y1 + y2);
            C(su1+1:end, :) = mchop(y3 + y4);
        else
            C = mchop(mchop(A.D) * mchop(B));
        end

    elseif isa(B, 'mphodlr')
        set_prec(prec_settings{B.level});

        if isempty(B.D)
            [mB, nB, su1, su2, sv1, sv2] = hsize(B);
            [mA, nA] = size(A);
            
            if nA ~= mB
                error('Please enter the inputs with correct dimensions.')
            end
        
            C = zeros(mA, nB);
            y1 = mphdot_dense(A(:, 1:su1), B.A11, prec_settings);
            y2 = mchop(mchop(mchop(A(:, su1+1:end)) * mchop(B.U2)) * mchop(B.V1));
            y3 = mchop(mchop(mchop(A(:, 1:su1)) * B.U1) * B.V2);
            y4 = mphdot_dense(A(:, su1+1:end), B.A22, prec_settings);
            
            C(:, 1:sv1) = mchop(y1 + y2);
            C(:, sv1+1:end) = mchop(y3 + y4);
        else
            C = mchop(mchop(A) * mchop(B.D));
        end
        
    else
        C = mchop(mchop(A) * mchop(B));
    end
end