function [t, emax] = base(prec)
    %Return the precision and maximum exponent
    if ismember(prec, {'q43','fp8-e4m3'})
        t = 4;
        emax = 7;
    elseif ismember(prec, {'q52','fp8-e5m2'})
        t = 3;
        emax = 15;
    elseif ismember(prec, {'h','half','fp16'})
        t = 11;
        emax = 15;
    elseif ismember(prec, {'b','bfloat16'})
        t = 8;
        emax = 127;  
    elseif ismember(prec, {'s','single','fp32'})
        t = 24;
        emax = 127;
    elseif ismember(prec, {'d','double','fp64'})
        t = 53;
        emax = 1023;
    else
        error('Please enter valid prec value.');
    end

end