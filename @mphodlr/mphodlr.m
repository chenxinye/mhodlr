classdef mphodlr
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
        The method to perform compression for off-diagonal blocks.

    threshold - double, default=1.0e-12
        The threshold value used for truncation of low rank approximation.

    type - str, default='dense'
        Under developed, used for detemine the HODLR matrix type.
    
    
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
        type {mustBeText} = 'dense' % TO DO
        prec_settings
        shape {mustBeNumeric} = []
        max_level {mustBeInteger} = 20
        bottom_level {mustBeInteger} = 0
    end

    properties(Access=private)
        min_block_size {mustBeInteger} = 20
        method {mustBeText} = 'svd'
        threshold  {mustBeNonNan, mustBeFinite, mustBeNumeric} = 1.0e-12
    end

    methods(Access=public)
        function obj = mphodlr(precs, A, varargin)
            obj.prec_settings = precs;

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
                obj.threshold = varargin{4};
            
            elseif nargin == 7
                obj.max_level = varargin{1};
                obj.min_block_size = varargin{2};
                obj.method = varargin{3};
                obj.threshold = varargin{4};
                obj.type = varargin{5};

            elseif nargin > 7
                disp(['Please enter the correct number or type of' ...
                    ' parameters.']);
            end
            
            obj.level = 1;
            min_size = min(size(A));
            [~, exponent] = log2(abs(min_size));
            
            if exponent < obj.max_level + 1
                obj.max_level = exponent - 1;
            end
            
            obj.check_exception();
            obj = build_hodlr_mat(obj, A, obj.level);
        end
        
        function obj =  build_hodlr_mat(obj, A, level)
            [rowSize, colSize] = size(A);
            
            obj.shape(1) = rowSize; 
            obj.shape(2) = colSize;
            
            if rowSize <= obj.min_block_size | colSize <= obj.min_block_size | level > obj.max_level
                obj.D = A;
                obj.bottom_level = max(obj.bottom_level, level-1);
                return;
            else
                obj.level = level;
                
                level = level + 1;
                rowSplit = ceil(rowSize / 2);
                colSplit = ceil(colSize / 2);

                obj.A11 = build_hodlr_mat(obj, A(1:rowSplit, 1:colSplit), ...
                    level);
                obj.A22 = build_hodlr_mat(obj, A(rowSplit+1:end, colSplit+1:end), ...
                    level);
                
                obj.bottom_level = max(obj.A11.bottom_level, obj.A22.bottom_level);
                if obj.level <= length(obj.prec_settings)
                    set_prec(obj.prec_settings{obj.level});
                    [obj.U1, obj.V2] = mp_compress(obj, A(1:rowSplit, colSplit+1:end));
                    [obj.U2, obj.V1] = mp_compress(obj, A(rowSplit+1:end, 1:colSplit));
                else
                    [obj.U1, obj.V2] = compress(obj, A(1:rowSplit, colSplit+1:end));
                    [obj.U2, obj.V1] = compress(obj, A(rowSplit+1:end, 1:colSplit));
                end
            end
        end
        

        function obj = transpose(obj)
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
                        'Approximation threshold: %d\n', obj.threshold);
                end
            end

            md = obj.method;
            td = obj.threshold;
            mbs = obj.min_block_size;
            ml = obj.max_level;
            tp = obj.type;

            if nargout > 1
                varargout = {md, td, mbs, ml, tp};
            end
        end
    end
   
    methods(Access=private)
        function [U, V] = compress(obj, A)
            [U, V] = compress_m(A, obj.method, obj.threshold);
        end

        function [U, V] = mp_compress(obj, A)
            [U, V] = mp_compress_m(A, obj.method, obj.threshold);
        end
        
        function check_exception(obj)

            if length(obj.prec_settings) < obj.max_level - 1
                warning(['The number of precisions used are less than ' ...
                    'the maximum tree level that can achieve. The remaining' ...
                    ' level will use the working precision for compresion.']); 
            end 
        end
    end

end

       

        