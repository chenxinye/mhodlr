function y = hstorage(obj, varargin)
    y = 0;
    if nargin == 1
        level = 1;
    else
        level = varargin{1};
    end

    if isa(obj, 'hodlr')
        if ~isempty(obj.D)
            if isa(obj.D, 'double') 
                bits = 64;
            else
                bits = 32;
            end

            [m, n] = size(obj.D);
            y = (m * n) * bits;
            return;
        else
            if isa(obj.U1, 'double') 
                bits = 64;
            else
                bits = 32;
            end

            [m, k] = size(obj.U1);
            [n, k] = size(obj.V2);
            
            y = y + (m + n) * k * bits;

            [m, k] = size(obj.U2);
            [n, k] = size(obj.V1);
            
            y = y + (m + n) * k * bits;
            y = y + hstorage(obj.A11, level+1);
            y = y + hstorage(obj.A22, level+1);
        end
    end
    
    if isa(obj, 'mphodlr')
        if ~isempty(obj.D)
            if isa(obj.D, 'double') 
                bits = 64;
            else
                bits = 32;
            end

            [m, n] = size(obj.D);
            y = (m * n) * bits;
            return;
        else
            if size(obj.prec_settings, 2) >= level
                bits = obj.prec_settings{level}.bits;
            else 
                bits = 64;
            end

            [m, k] = size(obj.U1);
            [n, k] = size(obj.V2);
            
            y = y + (m + n) * k * bits;

            [m, k] = size(obj.U2);
            [n, k] = size(obj.V1);
            
            y = y + (m + n) * k * bits;
            y = y + hstorage(obj.A11, level+1);
            y = y + hstorage(obj.A22, level+1);
        end
    end

    if isa(obj, 'amphodlr')
        if ~isempty(obj.D)
            if isa(obj.D, 'double') 
                bits = 64;
            else
                bits = 32;
            end

            [m, n] = size(obj.D);
            y = (m * n) * bits;
            return;
        else
            if obj.precIndex(level) ~= 0
                bits = obj.prec_settings{obj.precIndex(level)}.bits;
            else 
                bits = 64;
            end

            [m, k] = size(obj.U1);
            [n, k] = size(obj.V2);
            
            y = y + (m + n) * k * bits;

            [m, k] = size(obj.U2);
            [n, k] = size(obj.V1);
            
            y = y + (m + n) * k * bits;
            y = y + hstorage(obj.A11, level+1);
            y = y + hstorage(obj.A22, level+1);
        end
    end
end