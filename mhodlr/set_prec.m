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
    opt.built_in = prec.built_in;
    opt.ftp = prec.ftp;
end
