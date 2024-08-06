function x = mchop(x)
    global opt;

    if isempty(opt)
        precision();
    end

    emin = 1-opt.emax;
    xmin = 2^emin; 
    emins = emin + 1 - opt.t; 
    xmins = 2^emins;
    xmax = 2^opt.emax * (2-2^(1-opt.t));

    [~, exponent] = log2(abs(x));
    exponent = exponent - 1;
    ktemp = (exponent < emin & exponent >= emins);
    
    if opt.explim
        k_sub = find(ktemp);
        k_norm = find(~ktemp);
    else 
        k_sub = [];
        k_norm = 1:length(ktemp(:));
    end   

    x(k_norm) = pow2(roundit(pow2(x(k_norm), opt.t-1-exponent(k_norm))), exponent(k_norm)-(opt.t-1));

    if ~isempty(k_sub)
        t1 = opt.t - max(emin-exponent(k_sub),0);
        x(k_sub) = pow2(roundit(pow2(x(k_sub), t1-1-exponent(k_sub))), exponent(k_sub)-(t1-1));
     end

    if opt.explim
        switch(opt.round)
            case {1,6}
                xboundary = 2^opt.emax * (2-(1/2)*2^(1-opt.t));
                x(find(x >= xboundary)) = inf;   
                x(find(x <= -xboundary)) = -inf; 

            case 2
                x(find(x > xmax)) = inf;
                x(find(x < -xmax & x ~= -inf)) = -xmax;

            case 3
                x(find(x > xmax & x ~= inf)) = xmax;
                x(find(x < -xmax)) = -inf;

            case {4,5}
                x(find(x > xmax & x ~= inf)) = xmax;
                x(find(x < -xmax & x ~= -inf)) = -xmax;
            end

        if opt.subnormal == 0
            min_rep = xmin;
        else
            min_rep = xmins;
        end

        k_small = abs(x) < min_rep;

        switch(opt.round)
            case 1
                if opt.subnormal == 0
                    k_round = k_small & abs(x) >= min_rep/2;
                else
                    k_round = k_small & abs(x) > min_rep/2;
                end

                x(k_round) = sign(x(k_round)) * min_rep;
                x(k_small & ~k_round) = 0;
            
            case 2
                k_round = k_small & x > 0 & x < min_rep;
                x(k_round) = min_rep;
                x(k_small & ~k_round) = 0;
            
            case 3
                k_round = k_small & x < 0 & x > -min_rep;
                x(k_round) = -min_rep;
                x(k_small & ~k_round) = 0;

            case {4,5,6}
                x(k_small) = 0;
        end
    end
end


function y = roundit(x)
    global opt;

    default_sign = @(x) sign(x) + (x==0); 
    
    switch opt.round
      
      case 1
        y = abs(x);
        u = round(y - (rem(y,2) == 0.5));
        u(find(u == -1)) = 0;
        y = sign(x).*u; 
    
      case 2
        y = ceil(x); 
    
      case 3
        y = floor(x); 
    
      case 4
        y = (x >= 0 | x == -inf) .* floor(x) + (x < 0 | x == inf) .* ceil(x);
    
      case {5, 6}
    
        y = abs(x); 
        frac = y - floor(y);
        k = find(frac ~= 0);
        if isempty(k)
           y = x; 
        else   
          rnd = opt.randfunc(length(k));
          vals = frac(k);  vals = vals(:);
    
          switch opt.round
            case 5 
                   j = (rnd <= vals);
            case 6       
                   j = (rnd <= 0.5);
          end      
          y(k(j)) = ceil(y(k(j)));
          y(k(~j)) = floor(y(k(~j)));
          y = default_sign(x).*y; 
       end   
       
      otherwise
        error('Unsupported value of opt.round.')  
                   
    end
    
    if opt.flip
        
       temp = rand(size(y));
       k = find(temp <= opt.prob); 
       if ~isempty(k)
          u = abs(y(k));
          b = randi(opt.probarams(1)-1,size(u,1),size(u,2));
          u = bitxor(u,2.^(b-1));
          y(k) = default_sign(y(k)).*u; 
       end
    end
end