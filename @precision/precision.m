classdef precision
    properties
        t {mustBeInteger} = 11
        emax {mustBeInteger} = 15
        round {mustBeInteger} = 1
        subnormal = 1
        explim = 1
        prob {mustBeNonNan, mustBeFinite, mustBeNumeric} = 0.5
        flip = 0
        randfunc = @(n) rand(n, 1)
    end

    methods
        function obj = precision(varargin)
            switch nargin
                case 0 
                    [obj.t, obj.emax] = fpbase('h');
                    obj.round = 1;
                    obj.subnormal = 1;
                    obj.explim = 1;
                    obj.prob = 0.5;
                    obj.flip = 0;
                    obj.randfunc = @(n) rand(n, 1);
                    
                case 1
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = 1;
                    obj.subnormal = 1;
                    obj.explim = 1;
                    obj.prob = 0.5;
                    obj.flip = 0;
                    obj.randfunc = @(n) rand(n, 1);

                case 2
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = 1;
                    obj.explim = 1;
                    obj.prob = 0.5;
                    obj.flip = 0;
                    obj.randfunc = @(n) rand(n, 1);

                case 3
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = 1;
                    obj.prob = 0.5;
                    obj.flip = 0;
                    obj.randfunc = @(n) rand(n, 1);

                case 4
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = varargin{4};
                    obj.prob = 0.5;
                    obj.flip = 0;
                    obj.randfunc = @(n) rand(n, 1);


                case 5
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = varargin{4};
                    obj.prob = varargin{5};
                    obj.flip = 0;
                    obj.randfunc = @(n) rand(n, 1);
                
                case 6
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = varargin{4};
                    obj.prob = varargin{5};
                    obj.flip = varargin{6};
                    obj.randfunc = @(n) rand(n, 1);
                    
                case 7
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = varargin{4};
                    obj.prob = varargin{5};
                    obj.flip = varargin{6};
                    obj.randfunc = varargin{7};

                case 8
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = varargin{4};
                    obj.prob = varargin{5};
                    obj.flip = varargin{6};
                    obj.randfunc = varargin{7};
            end
        
            if nargin >= 1
                if strcmp(class(varargin{1}), 'double')
                    if length(varargin{1}) < 1 | length(varargin{1}) > 2
                        error('Please enter an valid value.');
                    else
                        obj.t = varargin{1};
                        obj.emax = 15;
                        return 
                    end
        
                    obj.t = varargin{1}(1);
                    obj.emax = varargin{1}(2);
        
                elseif ismember(varargin{1},  {'h','half','fp16','b','bfloat16','s', ...
                    'single','fp32','d','double','fp64',...
                    'q43','fp8-e4m3','q52','fp8-e5m2'})
                    [t, emax] = base(varargin{1});
                    obj.t = t;
                    obj.emax = emax;
                else
                    error('Please enter an valid value.');
                end
            end
        end
    end
end


