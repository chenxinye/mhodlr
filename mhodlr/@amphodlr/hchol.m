function R = hchol(H, varargin)
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
    if nargin < 2
        oformat = 'hodlr';
    else
        oformat = varargin{1};
    end
    

    if strcmp(oformat, 'dense')
        if isempty(H.D)
            [m, n, m1, m2, n1, n2] = hsize(H);
            R11 = hchol(H.A11, varargin{1});
            R12 = mldivide(R11', H.U1 * H.V2);
            R22 = hchol(fusedma(H.A22, -R12', R12, H.vareps), oformat);
            
            R = [R11, R12; zeros(m2, n1), R22];
        else
            R = chol(H.D);
        end
    else
        R = H;

        if isempty(H.D)
            R.A11 = hchol(H.A11, oformat);
            R12 = htrsl(R.A11.transpose(), H.U1 * H.V2);
            [R.U1, R.V2] = compress_m(R12, H.method, H.vareps, H.issparse);
            R.A22 = hchol(fusedma(H.A22, -R12', R12, H.vareps), oformat);
            R.U2 = zeros(size(R.U2,1),0);
            R.V1 = zeros(0,size(R.V1,2));
        else
            R.D = chol(H.D);
        end

    end
end


