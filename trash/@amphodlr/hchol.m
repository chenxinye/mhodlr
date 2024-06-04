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
        otype = 0;
    else
        if strcmp(varargin{1}, 'dense')
            otype = 1;
        else
            otype = 0;
        end    
    end
    

    if otype == 1
        if isempty(H.D)
            [m, n, m1, m2, n1, n2] = hsize(H);
            R11 = hchol(H.A11, H.threshold);
            R12 = mldivide(R11', H.U1 * H.V2);
            R22 = hchol(hrank_update(H.A22, -R12', R12, H.threshold), H.threshold);
            
            R = [R11, R12; zeros(m2, n1), R22];
        else
            R = chol(H.D);
        end
    else
        R = H;

        if isempty(H.D)
            R.A11 = hchol(H.A11, H.threshold);
            R12 = htrsl(R.A11.transpose(), H.U1 * H.V2);
            [R.U1, R.V2] = compress_m(R12, H.method, H.threshold);
            R.A22 = hchol(hrank_update(H.A22, -R12', R12, H.threshold), H.threshold);
            R.U2 = zeros(size(R.U2,1),0);
            R.V1 = zeros(0,size(R.V1,2));
        else
            R.D = chol(H.D);
        end

    end
end

%{  

    Compare below

function R = hchol(H, varargin)


    if nargin < 2
        epsilon = 1.0e-12;

    if isempty(H.D)
        [m, n, m1, m2, n1, n2] = hsize(H);
        R11 = hchol(H.A11);
        R12 = mldivide(R11', H.U1*H.V2);
        R22 = hchol(hadd_partial_hodlr(H.A22, R12' * R12, '-'));
        R = [R11, R12; zeros(m2, n1), R22];
    else
        R = chol(H.D);
    end
end



    %}

