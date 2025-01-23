function R = mhchol(H, varargin)
%{
    Compute Cholesky factorization for symmetric positive-definite HODLR matrix H.

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
    
    oformat - str, default='hodlr'
        The output format, either 'hodlr' or ''dense.


    Returns
    --------------------
    R - double
        The upper triangular matrix R is computed such that R' * R = H.  
%}
    if nargin < 3
        oformat = 'hodlr';
    else
        oformat = varargin{1};
    end
    
    if H.level <= 1
        global opt;
        set_prec(opt);
    end 

    H = hmchop(H);

    if strcmp(oformat, 'dense')
        if isempty(H.D)
            [m, n, m1, m2, n1, n2] = size_t(H);
            R11 = mhchol(H.A11, oformat);
            R12 = mchop(mldivide(R11', mchop(H.U1 * H.V2)));
            R22 = mhchol(mfusedma(H.A22, -R12', R12, H.vareps), oformat);
            
            R = [R11, R12; zeros(m2, n1), R22];
        else
            R = mchop(chol(H.D));
        end
    else
        R = H;

        if isempty(H.D)
            R.A11 = mhchol(H.A11, oformat);
            R12 = mhtrsl(R.A11.transpose(), mchop(H.U1 * H.V2));
            [R.U1, R.V2] = compress_m(R12, H.method, H.vareps, H.max_rnk, H.trun_norm_tp, H.issparse);
            R.A22 = mhchol(fusedma(H.A22, -R12', R12, H.vareps), oformat);
            R.U2 = zeros(size(R.U2,1),0);
            R.V1 = zeros(0,size(R.V1,2));
        else
            R.D = mchop(chol(H.D));
        end

    end
end


