function C = hadd(varargin)

    if nargin == 2
        operator = '+';
        oformat = 'hodlr';

    elseif nargin == 3
        operator = varargin{3};
        oformat = 'hodlr';

    elseif nargin == 4 
        operator = varargin{3};
        oformat = varargin{4};

    elseif nargin > 4
        error('Please enter the correct number of inputs.');
    end
    
    [mb, nb] = size(varargin{2});
    

    if strcmp(oformat, 'hodlr')
        C = varargin{1};

        if isempty(varargin{1}.D) 
            [m, n, m1, m2, n1, n2] = hsize(varargin{1});

            if m ~= mb | n ~= nb
                error('Please enter the inputs with consistent dimensions.');
            end  
            
            C.A11 = hadd(varargin{1}.A11, varargin{2}(1:m1, 1:n1), operator);
            C.A22 = hadd(varargin{1}.A22, varargin{2}(m1+1:end, n1+1:end), operator);

            if operator == '+'
                [C.U2, C.V1] = varargin{1}.compress(varargin{1}.U2 * varargin{1}.V1 + varargin{2}(m1+1:end, 1:n1));
                [C.U1, C.V2] = varargin{1}.compress(varargin{1}.U1 * varargin{1}.V2 + varargin{2}(1:m1, n1+1:end));
            else
                [C.U2, C.V1] = varargin{1}.compress(varargin{1}.U2 * varargin{1}.V1 - varargin{2}(m1+1:end, 1:n1));
                [C.U1, C.V2] = varargin{1}.compress(varargin{1}.U1 * varargin{1}.V2 - varargin{2}(1:m1, n1+1:end));
            end
        else
            if operator == '+'
                C.D = varargin{1}.D + varargin{2};
            else
                C.D = varargin{1}.D - varargin{2};
            end
        end
    else
        C = zeros(mb, nb);
        
        if isempty(varargin{1}.D) 
            [m, n, m1, m2, n1, n2] = hsize(varargin{1});
            
            if m ~= mb | n ~= nb
                error('Please enter the inputs with consistent dimensions.');
            end  

            C(1:m1, 1:n1) = hadd(varargin{1}.A11, varargin{2}(1:m1, 1:n1), operator, oformat);
            C(m1+1:end, n1+1:end) = hadd(varargin{1}.A22, varargin{2}(m1+1:end, n1+1:end), operator, oformat);
            
            if operator == '+'
                C(m1+1:end, 1:n1) = varargin{1}.U2 * varargin{1}.V1 + varargin{2}(m1+1:end, 1:n1);
                C(1:m1, n1+1:end) = varargin{1}.U1 * varargin{1}.V2 + varargin{2}(1:m1, n1+1:end);
            else
                C(m1+1:end, 1:n1) = varargin{1}.U2 * varargin{1}.V1 - varargin{2}(m1+1:end, 1:n1);
                C(1:m1, n1+1:end) = varargin{1}.U1 * varargin{1}.V2 - varargin{2}(1:m1, n1+1:end);
            end

        else
            if operator == '+'
                C = varargin{1}.D + varargin{2};
            else
                C = varargin{1}.D - varargin{2};
            end
        end 
    end
end