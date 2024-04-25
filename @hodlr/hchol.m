function R = hchol(H, varargin)
    if nargin < 2
        epsilon = 1.0e-12;

    if isempty(H.D)
        [m, n, m1, m2, n1, n2] = hsize(H);
        R11 = hchol(H.A11);
        R12 = mldivide(R11', H.U1);
        X = -H.U2 * (R12' * R12);
        R22 = hchol(hrank_update(H.A22, X, H.V2, epsilon));
        
        R12 = R12 * H.V2;
        R = [R11, R12; zeros(m2, n1), R22];
    else
        R = chol(H.D);
    end
end

% Compare to below
% function R = hchol(H, varargin)
%     if nargin < 2
%         epsilon = 1.0e-12;
%
%     if isempty(H.D)
%         [m, n, m1, m2, n1, n2] = hsize(H);
%         R11 = hchol(H.A11);
%         R12 = mldivide(R11', H.U1*H.V2);
%        R22 = hchol(H.A22);
%        R = [R11, R12; zeros(m2, n1), R22];
%    else
%        R = chol(H.D);
%    end
%end