classdef hodlr
    properties(Access=public)
        A11
        A22
        
        U1
        V2

        U2
        V1
           
        D
        level
    end

    properties(Access=private)
        min_block_size
        max_level
        method
        threshold
    end

    methods(Access=public)
        function obj = hodlr(varargin)
            if nargin == 1
                obj.method = 'svd';
                obj.threshold = 1.0e-12;
                obj.min_block_size = 2;
                obj.max_level = 9999;

            elseif nargin == 2
                obj.method = varargin{2};
                obj.threshold = 1.0e-12;
                obj.min_block_size = 2;
                obj.max_level = 9999;

            elseif nargin == 3
                obj.method = varargin{2};
                obj.threshold = varargin{3};
                obj.min_block_size = 2;
                obj.max_level = 9999;

            elseif nargin == 4
                obj.method = varargin{2};
                obj.threshold = varargin{3};
                obj.min_block_size = varargin{4};
                obj.max_level = 9999;

            elseif nargin == 5
                obj.method = varargin{2};
                obj.threshold = varargin{3};
                obj.min_block_size = varargin{4};
                obj.max_level = varargin{5};

            else 
                disp(['Please enter the correct number or type of' ...
                    ' parameters.']);
            end
            
            obj.level = 1;
            obj = build_hodlr_mat(obj, varargin{1}, obj.level);
        end
        
        function obj =  build_hodlr_mat(obj, A, level)
            [rowSize, colSize] = size(A);
            
            if rowSize <= obj.min_block_size | colSize <= obj.min_block_size
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
                    
                [obj.U1, obj.V2] = compress(obj, A(1:rowSplit, colSplit+1:end));
                [obj.U2, obj.V1] = compress(obj, A(rowSplit+1:end, 1:colSplit));
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


        function C = inverse_hodlr(obj)
            if ~issquare(obj)
                error('Inverse is only applied to a square HODLR matrix.');
            end
            
            C = obj;
            if isempty(obj.D)
                X22 = inverse_hodlr(obj.A22);
                A12 = obj.U1 * obj.V2;
                A21 = obj.U2 * obj.V1;
                
                C.A11  = inverse_hodlr(hadd(obj.A11, hdot_double(hdot(A12, X22), A21), '-'));
            
                [C.U1, C.V2] = compress_m(hdot_double(hdot_double(C.A11, -A12), X22), obj.method, obj.threshold);
                C21 = -hdot_double(hdot_double(X22, A21), C.A11);
                [C.U2, C.V1] = compress_m(C21, obj.method, obj.threshold);
                XX = hdot_double(C21 * A12, X22);
                C.A22 = hadd(X22, XX, '-');
            else
                C.D = inv(obj.D);
            end
        end
        
        function C = inverse_double(obj)
            if ~issquare(obj)
                error('Inverse is only applied to a square HODLR matrix.');
            end
            
            if class(obj) == 'hodlr'
                if isempty(obj.D)
                    X22 = inverse_double(obj.A22);
                    A12 = obj.U1*obj.V2;
                    A21 = obj.U2*obj.V1;
                    X11 = inv(hadd_partial_double(obj.A11, A12 * X22 * A21 ,'-'));
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

        function load_params(obj)
            fprintf(...
                'minimum block size: %d\n', obj.min_block_size);
            fprintf(...
                'Approximation method: %s\n', obj.method);
            fprintf(...
                'Approximation threshold: %d\n', obj.threshold);
        end
    end

    methods(Access=private)
        function [U, V] = compress(obj, A)
            [U, V] = compress_m(A, obj.method, obj.threshold);
        end
    end

end

       

        