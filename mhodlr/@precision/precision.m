classdef precision
%{
    Parameters
    --------------------
    base - array | str, default='h'
       For the string type, the arithmetic format supports:
       'q43', 'fp8-e4m3'       - NVIDIA quarter precision (4 exponent bits,
                                 3 significand (mantissa) bits),
       'q52', 'fp8-e5m2'       - NVIDIA quarter precision (5 exponent bits,
                                 2 significand bits),
       'b', 'bfloat16'         - bfloat16,
       'h', 'half', 'fp16'     - IEEE half precision (the default),
       's', 'single', 'fp32'   - IEEE single precision,
       'd', 'double', 'fp64'   - IEEE double precision,
       'c', 'custom'           - custom format.

      The custom (base 2) format is defined by array (2-element, [t, emax]), 
      where t is the number of bits in the significand
      (including the hidden bit) and emax is the maximum value of the
      exponent.  The minimu exponent is taken to be emin = 1 - emax and
      the IEEE floating-point number representation is assumed, so that
      emax and the number of  bits e in the exponent are related by
      emax = 2^(e-1) - 1. 

    round - int, default=1
       1: round to nearest using round to even last bit to break ties
          (the default);
       2: round towards plus infinity (round up);
       3: round towards minus infinity (round down);
       4: round towards zero;
       5: stochastic rounding - round to the next larger or next smaller
          f.p. (floating-point) number with probability proportional to
          1 minus the distance to those f.p. numbers;
       6: stochastic rounding - round to the next larger or next smaller 
          f.p. number with equal probability.

    subnormal - int
        specifies whether subnormal numbers are supported, if `subnormal=0`, 
        subnormals are flushed to zero:
            0 = do not support subnormals (the default for base='b', i.e., bfloat16 format),
            1 = support subnormals (the default if base is set to others).

    explim - int, default=1
        ``explim = 0`` make emax (the maximal exponent) for the specified arithmetic disabled, 
        so overflow, underflow, or subnormal numbers will be produced only if necessary 
        for the data type.  This parameter is for exploring
        low precisions independent of range limitations.

    flip - int, default=0
        Determine whether each element of the rounded value has a randomly chosen bit in 
        its significand flipped with a certain probability.

    prob - double, default=0.5
         if flip = 1 then each element of the rounded
        value has a randomly chosen bit in its significand flipped with probability ``prob``.

    randfunc - func, default=@(n) rand(n, 1)
        The random function for stochastic rounding. 
        If options.randfunc is supplied, then in stochastic rounding (modes
        5 and 6) the random numbers used for rounding will be generated
        using that function. It should be a function that has a single argument
        for the number of random numbers to generate and returns a vector of
        the random numbers. 
        

    Properties
    --------------------
    Same as parameters, except 

    u - double
        Unit roundoff computed for current floating point format.
%}

    properties
        t {mustBeInteger} = 11
        emax {mustBeInteger} = 15
        round {mustBeInteger} = 1
        subnormal = 1
        explim = 1
        prob {mustBeNonNan, mustBeFinite, mustBeNumeric} = 0.5
        flip = 0
        randfunc = @(n) rand(n, 1)
        u {mustBeFinite, mustBeNumeric}
        bits
        ftp = 'custom'
        built_in = false
    end

    methods
        function obj = precision(varargin)
            switch nargin
                case 0 
                    [obj.t, obj.emax] = base('h');
                    obj.round = 1;
                    obj.subnormal = 1;
                    obj.explim = 1;
                    obj.flip = 0;
                    obj.prob = 0.5;
                    obj.randfunc = @(n) rand(n, 1);
                    
                case 1
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = 1;
                    obj.subnormal = 1;
                    obj.explim = 1;
                    obj.flip = 0;
                    obj.prob = 0.5;
                    obj.randfunc = @(n) rand(n, 1);

                case 2
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = 1;
                    obj.explim = 1;
                    obj.flip = 0;
                    obj.prob = 0.5;
                    obj.randfunc = @(n) rand(n, 1);

                case 3
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = 1;
                    obj.flip = 0;
                    obj.prob = 0.5;
                    obj.randfunc = @(n) rand(n, 1);

                case 4
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = varargin{4};
                    obj.flip = 0;
                    obj.prob = 0.5;
                    obj.randfunc = @(n) rand(n, 1);


                case 5
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = varargin{4};
                    obj.flip = varargin{5};
                    obj.prob = 0.5;
                    obj.randfunc = @(n) rand(n, 1);
                
                case 6
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = varargin{4};
                    obj.flip = varargin{5};
                    obj.prob = varargin{6};
                    obj.randfunc = @(n) rand(n, 1);
                    
                case 7
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = varargin{4};
                    obj.flip = varargin{5};
                    obj.prob = varargin{6};
                    obj.randfunc = varargin{7};
                
                otherwise
                    obj.t = 0;
                    obj.emax = 0;
                    obj.round = varargin{2};
                    obj.subnormal = varargin{3};
                    obj.explim = varargin{4};
                    obj.flip = varargin{5};
                    obj.prob = varargin{6};
                    obj.randfunc = varargin{7};
            end

            if nargin >= 1
                if strcmp(class(varargin{1}), 'double')
                    if length(varargin{1}) < 1 | length(varargin{1}) > 2
                        error('Please enter an valid value.');
                    elseif length(varargin{1}) == 1
                        obj.t = varargin{1};
                        obj.emax = 15;
                        obj.u = 2^(1 - obj.t) / 2; 
                        obj.bits = obj.t + log2((obj.emax + 1)* 2);
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
                    
                    obj.ftp = varargin{1};

                    if ismember(varargin{1}, {'b','bfloat16'})
                        obj.subnormal = 0;
                    end
                else
                    error('Please enter an valid value.');
                end
            end

            obj.u = 2^(1 - obj.t) / 2; 
            obj.bits = obj.t + round(log2((obj.emax + 1)* 2));
        end


        function obj = builtin(obj, varargin)
            if ~ismember(obj.ftp, {'h','half','fp16','s', 'single','fp32','d','double','fp64'})
                error("The current floating point types does not support built-in rounding.")
            end

            if nargin == 1
                obj.built_in = true;
            else
                obj.built_in = varargin{1};
            end
        end
    end
end


