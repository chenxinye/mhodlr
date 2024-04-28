function R = hchol(H, varargin)
%% Compute Cholesky factorization for symmetric positive-definite HODLR matrix H.
%
% Parameters
% --------------------
% H - hodlr
%     Matrix in HODLR format - hodlr class.
% 
% epsilon - double
%     The threshold for recompression.
%
%
% Returns
% --------------------
% R - double
%     The upper triangular matrix R is computed such that R' * R = H. 
% 
    if nargin < 2
        epsilon = 1.0e-12;

    if isempty(H.D)
        [m, n, m1, m2, n1, n2] = hsize(H);
        R11 = hchol(H.A11);
        R12 = mldivide(R11', H.U1 * H.V2);
        R22 = hchol(hrank_update(H.A22, -R12', R12, epsilon));
        
        R = [R11, R12; zeros(m2, n1), R22];
    else
        R = chol(H.D);
    end
end

%{  

% Compare below

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

