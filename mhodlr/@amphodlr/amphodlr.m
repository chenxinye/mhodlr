classdef amphodlr
%{
    Parameters
    --------------------
    precs - cell
        The cell array that contains the precision used for compression of each level. 
        Each element is a precision class.

    A - double | single
        Matrix to be converted.
        
    max_level - int, default=9999
        Maximum level for cluster tree.

    min_block_size - int, default=2
        The minimum size for HODLR blocks.

    method - str, default='svd'
        The method to perform compression for off-diagonal blocks. Options: 'rsvd' and 'qr' used for randomized SVD or QR factorization.

    vareps - double, default=1.0e-12
        The vareps value used for truncation of low rank approximation.
    
    max_rnk - int, default=999
        The maximum rank of the off-diagonal block used for low rank truncation.

    trun_norm_tp - str, default='2'
        Norm type for the the off-diagonal block truncation ``||A - B||_trun_norm_tp <= vareps * ||B||``.
        
    issparse - bool, default=false:
        Whether or not store the generators U and V in sparse format.
    
    Properties
    --------------------
    U1, V2, U2, V1 - double 
        The right upper block matrix of each level, we have A12 = U1 * V2 and A21 = U2 * V1.

    A11, A22 - hodlr 
        The diagonal block matrix in HODLR format (access in the next level). 

    shape - array
        The shape of object in the current level. 

    level - int
        The level for cluster tree.
    
    max_level - int
        The maximum level of cluster tree after transformation.

%}

    properties(Access=public)

        U1 {mustBeNonNan, mustBeFinite, mustBeNumeric}
        V2 {mustBeNonNan, mustBeFinite, mustBeNumeric}

        U2 {mustBeNonNan, mustBeFinite, mustBeNumeric}
        V1 {mustBeNonNan, mustBeFinite, mustBeNumeric}
           
        D {mustBeNonNan, mustBeFinite, mustBeNumeric}

        A11 % HODLR format
        A22 % HODLR format
        
        level {mustBeInteger} = 0
        
        shape {mustBeNumeric} = []
        max_level {mustBeInteger} = 20
        bottom_level {mustBeInteger} = 0
        normOrder {mustBeNonNan, mustBeFinite, mustBeNumeric}
        precIndex {mustBeNonNan, mustBeFinite, mustBeNumeric}
        unitRoundOff {mustBeNonNan, mustBeFinite, mustBeNumeric}

        min_block_size {mustBeInteger} = 20
        vareps {mustBeNonNan, mustBeFinite, mustBeNumeric} = 1.0e-12
        max_rnk = 999
        prec_settings
        issparse = false
    end

    properties(Access=private)
        method {mustBeText} = 'svd'
        precIndexBool {mustBeNonNan, mustBeFinite}
        sortIdx
        trun_norm_tp = '2'
    end

    methods(Access=public)
        function obj = amphodlr(precs, A, varargin)
            if nargin == 0
                obj.D = [];
                
            elseif strcmp(class(precs), 'char')
                if strcmp(precs, 'eye')
                    A = eye(A);
                    obj.max_level = varargin{1};
          
                    if nargin == 4
                        obj.min_block_size = varargin{2};
                    end

                    obj = build_hodlr_eye(obj, A, 1);

                elseif strcmp(precs, 'ones')
                    A = ones(A);
                    obj.max_level = varargin{1};

                    if nargin == 4
                        obj.min_block_size = varargin{2};
                    end
                    
                    obj = build_hodlr_ones(obj, A, 1);
                
                elseif strcmp(precs, 'zeros')
                    A = zeros(A);
                    obj.max_level = varargin{1};

                    if nargin == 4
                        obj.min_block_size = varargin{2};
                    end

                    obj = build_hodlr_zeros(obj, A, 1);
                end
                  
            else
                obj.prec_settings = [prec_chain(precision('d')), precs];
                
                [obj.unitRoundOff, obj.sortIdx] = sort_by_u(obj, obj.prec_settings);
                % obj.prec_settings = obj.prec_settings(obj.sortIdx);
                
                if nargin == 3
                    obj.max_level = varargin{1};

                elseif nargin == 4
                    obj.max_level = varargin{1};
                    obj.min_block_size = varargin{2};

                elseif nargin == 5
                    obj.max_level = varargin{1};
                    obj.min_block_size = varargin{2};
                    obj.method = varargin{3};

                elseif nargin == 6
                    obj.max_level = varargin{1};
                    obj.min_block_size = varargin{2};
                    obj.method = varargin{3};
                    obj.vareps = varargin{4};
                
                elseif nargin == 7
                    obj.max_level = varargin{1};
                    obj.min_block_size = varargin{2};
                    obj.method = varargin{3};
                    obj.vareps = varargin{4};
                    obj.max_rnk = varargin{5};

                elseif nargin == 8
                    obj.max_level = varargin{1};
                    obj.min_block_size = varargin{2};
                    obj.method = varargin{3};
                    obj.vareps = varargin{4};
                    obj.max_rnk = varargin{5};
                    obj.trun_norm_tp = varargin{6};

                elseif nargin == 9
                    obj.max_level = varargin{1};
                    obj.min_block_size = varargin{2};
                    obj.method = varargin{3};
                    obj.vareps = varargin{4};
                    obj.max_rnk = varargin{5};
                    obj.trun_norm_tp = varargin{6};
                    obj.issparse = varargin{6};

                elseif nargin > 8
                    disp(['Please enter the correct number or type of' ...
                        ' parameters.']);
                end
                
                if strcmp(class(A) , 'single')
                    obj.issparse = false;
                end

                obj.level = 1;
                min_size = min(size(A));
                max_level = floor(log2(abs(min_size)));

                if obj.max_level > max_level
                    obj.max_level = max_level;
                end

                obj.normOrder = zeros(1, obj.max_level+1);
                obj.precIndex = ones(1, obj.max_level);
                obj.precIndexBool = zeros(1, obj.max_level);
                obj.normOrder(1) = sum(A.^2, 'all');
                [obj, obj.normOrder] = initialize(obj, A, obj.level, obj.normOrder);
                [obj, obj.precIndex, obj.precIndexBool] = build_hodlr_mat(obj, A, obj.level, ...
                                                        obj.precIndex, obj.precIndexBool);

                obj.precIndex = obj.precIndex(1: obj.bottom_level);
                obj.precIndexBool = obj.precIndexBool(1: obj.bottom_level);
            end
        end


        function [obj, normOrder] =  initialize(obj, A, level, normOrder)
            [rowSize, colSize] = size(A);
            
            obj.shape(1) = rowSize;
            obj.shape(2) = colSize;
            
            if floor(rowSize / 2) < obj.min_block_size | floor(colSize / 2) < obj.min_block_size | level > obj.max_level
                return;
            else
                obj.level = level;

                level = level + 1;
                rowSplit = ceil(rowSize / 2);
                colSplit = ceil(colSize / 2);

                nrm1 = sum(A(1:rowSplit, colSplit+1:end).^2, 'all');
                nrm2 = sum(A(rowSplit+1:end, 1:colSplit).^2, 'all');
                
                normOrder(obj.level+1) = max([normOrder(obj.level+1), nrm1, nrm2]);
                
                [obj.A11, normOrder] = initialize(obj, A(1:rowSplit, 1:colSplit), level, normOrder);
                [obj.A22, normOrder] = initialize(obj, A(rowSplit+1:end, colSplit+1:end), level, normOrder);
                
            end
        end

        function obj =  build_hodlr_zeros(obj, A, level)
            [rowSize, colSize] = size(A);
            
            obj.shape(1) = rowSize; 
            obj.shape(2) = colSize;
           
            if floor(rowSize / 2) < obj.min_block_size | floor(colSize / 2) < obj.min_block_size | level > obj.max_level
                obj.D = A;
                obj.bottom_level = max(obj.bottom_level, level-1);
                return;
            else
                obj.level = level;
    
                level = level + 1;
                rowSplit = ceil(rowSize / 2);
                colSplit = ceil(colSize / 2);

                obj.A11 = build_hodlr_zeros(obj, A(1:rowSplit, 1:colSplit), ...
                    level);
                obj.A22 = build_hodlr_zeros(obj, A(rowSplit+1:end, colSplit+1:end), ...
                    level);
                
                obj.bottom_level = max(obj.A11.bottom_level, obj.A22.bottom_level);
                obj.U1 = zeros(rowSplit, 1);
                obj.V2 = zeros(1, colSize-colSplit);
                obj.U2 = zeros(rowSize-rowSplit, 1);
                obj.V1 = zeros(1, colSplit);
            end
        end

        function obj =  build_hodlr_eye(obj, A, level)
            [rowSize, colSize] = size(A);
            
            obj.shape(1) = rowSize; 
            obj.shape(2) = colSize;
            
            if floor(rowSize / 2) < obj.min_block_size | floor(colSize / 2) < obj.min_block_size | level > obj.max_level
                obj.D = A;
                obj.bottom_level = max(obj.bottom_level, level-1);
                return;
            else
                obj.level = level;
    
                level = level + 1;
                rowSplit = ceil(rowSize / 2);
                colSplit = ceil(colSize / 2);

                obj.A11 = build_hodlr_eye(obj, A(1:rowSplit, 1:colSplit), ...
                    level);
                obj.A22 = build_hodlr_eye(obj, A(rowSplit+1:end, colSplit+1:end), ...
                    level);
                
                obj.bottom_level = max(obj.A11.bottom_level, obj.A22.bottom_level);
                obj.U1 = zeros(rowSplit, 1);
                obj.V2 = zeros(1, colSize-colSplit);
                obj.U2 = zeros(rowSize-rowSplit, 1);
                obj.V1 = zeros(1, colSplit);
            end
        end

        function obj =  build_hodlr_ones(obj, A, level)
            [rowSize, colSize] = size(A);
            
            obj.shape(1) = rowSize; 
            obj.shape(2) = colSize;
           
            if floor(rowSize / 2) < obj.min_block_size | floor(colSize / 2) < obj.min_block_size | level > obj.max_level
                obj.D = A;
                obj.bottom_level = max(obj.bottom_level, level-1);
                return;
            else
                obj.level = level;
    
                level = level + 1;
                rowSplit = ceil(rowSize / 2);
                colSplit = ceil(colSize / 2);

                obj.A11 = build_hodlr_ones(obj, A(1:rowSplit, 1:colSplit), ...
                    level);
                obj.A22 = build_hodlr_ones(obj, A(rowSplit+1:end, colSplit+1:end), ...
                    level);
                
                obj.bottom_level = max(obj.A11.bottom_level, obj.A22.bottom_level);
                obj.U1 = ones(rowSplit, 1);
                obj.V2 = ones(1, colSize-colSplit);
                obj.U2 = ones(rowSize-rowSplit, 1);
                obj.V1 = ones(1, colSplit);
            end
        end


        function [obj, precIndex, precIndexBool] =  build_hodlr_mat(obj, A, level, ...
                                                        precIndex, precIndexBool)
            [rowSize, colSize] = size(A);
            
            obj.shape(1) = rowSize; 
            obj.shape(2) = colSize;
            
            if floor(rowSize / 2) < obj.min_block_size | floor(colSize / 2) < obj.min_block_size | level > obj.max_level
                obj.D = A;
                obj.bottom_level = max(obj.bottom_level, level-1);
                obj.max_rnk = min(rowSize, colSize);
                return;
            else
                obj.level = level;
                
                level = level + 1;
                rowSplit = ceil(rowSize / 2);
                colSplit = ceil(colSize / 2);

                if precIndexBool(obj.level) == 0
                    xi = sqrt(obj.normOrder(obj.level+1) / obj.normOrder(1)) ;
                    update_u = obj.vareps / (2^((obj.level+1)/2) * xi);
                    if ~isinf(update_u)
                        find_u = find(obj.unitRoundOff<=update_u);
    
                        if ~isempty(find_u)
                            % disp('---------------------')
                            % disp(obj.level)
                            precIndex(obj.level) = obj.sortIdx(find_u(end)) - 1;
                            % disp(precIndex(obj.level))
                        else
                            error('Not available precision');
                        end
                    else
                        precIndex(obj.level) = precIndex(obj.level-1);
                    end
                end

                precIndexBool(obj.level) = 1;

                [obj.A11, precIndex, precIndexBool] = build_hodlr_mat(obj, A(1:rowSplit, 1:colSplit), ...
                    level, precIndex, precIndexBool);
                
                [obj.A22, precIndex, precIndexBool] = build_hodlr_mat(obj, A(rowSplit+1:end, colSplit+1:end), ...
                    level, precIndex, precIndexBool);

                obj.bottom_level = max(obj.A11.bottom_level, obj.A22.bottom_level);

                set_prec(obj.prec_settings{precIndex(obj.level)+1});
                
                [obj.U1, obj.V2, max_rnk1] = mp_compress(obj, A(1:rowSplit, colSplit+1:end));
                [obj.U2, obj.V1, max_rnk2] = mp_compress(obj, A(rowSplit+1:end, 1:colSplit));
                
                if obj.level < obj.bottom_level - 1
                    obj.max_rnk = max([max_rnk1, max_rnk2]);
                else
                    obj.max_rnk = max([max_rnk1, max_rnk2, obj.A11.max_rnk, obj.A22.max_rnk]);
                end
            end
        end
        

        function obj = transpose(obj)
            temp = obj.shape(1);
            obj.shape(1) = obj.shape(2);
            obj.shape(2) = temp;
            
            if isempty(obj.D)
                copyU2 = obj.U2;
                copyV1 = obj.V1;
                obj.U2 = obj.V2.';
                obj.V1 = obj.U1.';

                obj.U1 = copyV1.';
                obj.V2 = copyU2.';
                obj.A11 = transpose(obj.A11);
                obj.A22 = transpose(obj.A22);
            else
                obj.D = obj.D.';
            end
        end

        function [A] = todense(obj)
            A = recover(obj);
        end

        function [varargout] = load_params(obj, varargin)
            %% Load parameters of HODLR matrix
            %% 
            % Parameters
            % --------------------
            % verbose - boolean, default=1
            %       Whether or not print the parameter settings. 

            if nargin == 2
                if varargin{1} == 1
                    fprintf(...
                        'Minimum block size: %d\n', obj.min_block_size);
                    fprintf(...
                        'Approximation method: %s\n', obj.method);
                    fprintf(...
                        'Tree depth: %d\n', obj.max_level);
                    fprintf(...
                        'Approximation vareps: %d\n', obj.vareps);
                end
            end

            md = obj.method;
            td = obj.vareps;
            mbs = obj.min_block_size;
            ml = obj.max_level;
            tp = obj.type;

            if nargout > 1
                varargout = {md, td, mbs, ml, tp};
            end
        end
    end
   
    methods(Access=private)
        function [U, V, max_rnk] = compress(obj, A)
            [U, V, max_rnk] = compress_m(A, obj.method, obj.vareps, obj.max_rnk, obj.trun_norm_tp, obj.issparse);
        end

        function [U, V, max_rnk] = mp_compress(obj, A)
            [U, V, max_rnk] = mp_compress_m(A, obj.method, obj.vareps, obj.max_rnk, obj.trun_norm_tp, obj.issparse);
        end
        
        function [sortu, sortIdx] = sort_by_u(obj, u_chain)
            callCellFunc = cellfun(@(x)x.u, u_chain);
            [sortu, sortIdx] = sort(callCellFunc);
        end
    end
end

       

        
