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
        The shape of input matrix A. 

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
        max_level {mustBeInteger} = 9999
    end

    properties(Access=private)
        min_block_size {mustBeInteger} = 2
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
            obj.shape(1) = size(A, 1); obj.shape(2) = size(A, 2);
            
            min_size = min(obj.shape);
            [~, exponent] = log2(abs(min_size));
            
            if exponent < obj.max_level + 1
                obj.max_level = exponent - 1;
            end
            
            obj.check_exception();
            obj = build_hodlr_mat(obj, A, obj.level);
        end
        
        function obj =  build_hodlr_mat(obj, A, level)
            [rowSize, colSize] = size(A);
            
            if rowSize <= obj.min_block_size | colSize <= obj.min_block_size | level >= obj.max_level
                obj.D = A;
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
                
                if obj.level < length(obj.prec_settings)
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

        function C = inverse(obj, varargin)
            %% The function is to return the inverse of HOLDR matrix.
            %
            % Parameters
            % --------------------
            % H - hodlr
            %     Matrix in HODLR format - hodlr class.
            % 
            % algorithm - int, default=1
            %     The algorithm to implement inverse
            % 
            % oformat - str, default = 'hodlr'
            %     The format of returns.
            % 
            % Returns
            % --------------------
            % C - hodlr | double
            %     Return matrix in hodlr class or dense array.
            % 
            
            switch nargin
                case 1 
                    algorithm = 1;
                    oformat = 'hodlr';

                case 2
                    algorithm = varargin{1};
                    oformat = 'hodlr';

                case 3 
                    algorithm = varargin{1};
                    oformat = varargin{2};
            end

            if algorithm == 1
                if algorithm == 'hodlr'
                    C = inverse_hodlr(obj);
                else
                    C = inverse_dense(obj);
                end
            else
                if algorithm == 'hodlr'
                    C = inverse_nonrecursive_hodlr(obj);
                else
                    C = inverse_nonrecursive_dense(obj);
                end
            end
        end 

        function C = inverse_nonrecursive_hodlr(obj)
            %% This member method implement inverse of hodlr matrix in nonrecursive manner
            [m, n] = hsize(obj);

            if m ~= n
                error('Inverse is only applied to a square HODLR matrix.');
            end

            if isempty(obj.D)
                C1 = inverse_dense(obj.A11);
                C2 = inverse_dense(obj.A22);
                
                A12 = obj.U1 * obj.V2;
                A21 = obj.U2 * obj.V1;
                X = blkoffdiag(C1 * A12, C2* A21);
                [U, D, V] = svd(X);

                L = eye(m) - U * inv(inv(D) + V' * U) * V';
                R = blkdiag(C1, C2);
                C = L * R;
            else
                C = inv(obj.D);
            end

            [md, td, mbs, ml, tp] = load_params(obj, 0);
            C = hodlr(C, md, td, mbs, ml, tp);
        end

        function C = inverse_nonrecursive_dense(obj)
            %% This member method implement inverse of hodlr matrix in nonrecursive manner
            [m, n] = hsize(obj);

            if m ~= n
                error('Inverse is only applied to a square HODLR matrix.');
            end
            
            if isempty(obj.D)
                C1 = inverse_nonrecursive_dense(obj.A11);
                C2 = inverse_nonrecursive_dense(obj.A22);
                
                A12 = obj.U1 * obj.V2;
                A21 = obj.U2 * obj.V1;
                X = blkoffdiag(C1 * A12, C2* A21);
                [U, D, V] = svd(X);

                L = eye(m) - U * inv(inv(D) + V' * U) * V';
                R = blkdiag(C1, C2);
                C = L * R;
            else
                C = inv(obj.D);
            end
            
        end

        function C = inverse_hodlr(obj)
            if ~issquare(obj)
                error('Inverse is only applied to a square HODLR matrix.');
            end
            
            C = obj;
            if isempty(obj.D)
                X22 = inverse_hodlr(obj.A22);
                A12 = obj.U1 * obj.V2;
                A21 = obj.U2 * obj.V1;
                
                C.A11  = inverse_hodlr(hadd(obj.A11, hdot_dense(hdot(A12, X22), A21), '-'));
            
                [C.U1, C.V2] = compress_m(hdot_dense(hdot_dense(C.A11, -A12), X22), obj.method, obj.threshold);
                C21 = -hdot_dense(hdot_dense(X22, A21), C.A11);
                [C.U2, C.V1] = compress_m(C21, obj.method, obj.threshold);
                XX = hdot_dense(C21 * A12, X22);
                C.A22 = hadd(X22, XX, '-');
            else
                C.D = inv(obj.D);
            end
        end
        
        function C = inverse_dense(obj)
            if ~issquare(obj)
                error('Inverse is only applied to a square HODLR matrix.');
            end
            
            if class(obj) == 'hodlr'
                if isempty(obj.D)
                    X22 = inverse_dense(obj.A22);
                    A12 = obj.U1*obj.V2;
                    A21 = obj.U2*obj.V1;
                    X11 = inverse_dense(hadd_partial_hodlr(obj.A11, A12 * X22 * A21 ,'-'));
                    C12 = -X11 * A12 * X22;
                    C21 = -X22 * A21 * X11;
                    C22 = X22 + X22 * A21 * X11 * A12 * X22;
    
                    C = [X11, C12; C21, C22];
                else
                    C = inv(obj.D);
                end
            else
                C = inv(obj);
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
                        'minimum block size: %d\n', obj.min_block_size);
                    fprintf(...
                        'Approximation method: %s\n', obj.method);
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

            if length(obj.prec_settings) < obj.max_level
                warning(['The number of precisions used are less than ' ...
                    'the maximum tree level that can achieve. The remaining' ...
                    ' level will use the working precision for compresion.']); 
            end 
        end
    end

end

       

        