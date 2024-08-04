function opt = set_prec(prec)

    global opt;

    opt.t = prec.t;
    opt.emax = prec.emax;
    opt.round = prec.round;
    opt.subnormal = prec.subnormal;
    opt.explim = prec.explim;
    opt.prob = prec.prob;
    opt.flip = prec.flip;
    opt.randfunc = prec.randfunc;
    
end