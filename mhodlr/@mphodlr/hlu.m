function [L, U] = hlu(H, varargin)
    %{
        Compute LU factorization for HODLR matrix H.
    
        Parameters
        --------------------
        H - hodlr
            Matrix in HODLR format - hodlr class.
        
        oformat - str, default='hodlr'
            The output format. 'dense' or 'hodlr'.
    
        epsilon - double, default is the vareps of holdlr matrix H
            The vareps for recompression.
    
        Returns
        --------------------
        L - double
            The upper triangular matrix L is computed such that L * U = H. 
        U - double
            The upper triangular matrix U is computed such that L * U = H. 
    %}
    
        if ~issquare(H)
            error('LU factorization is only applied to a square HODLR matrix.');
        end
        
        if nargin == 1
            oformat = 'hodlr';
            epsilon = H.vareps;
        elseif nargin == 2 
            oformat = varargin{1};
            epsilon = H.vareps;
        else
            oformat = varargin{1};
            epsilon = varargin{2};
        end 
        
        [m, n] = hsize(H);
        
    
        if strcmp(oformat, 'dense')
            if isempty(H.D)
                [L11, U11] = hlu(H.A11, 'dense', epsilon);
                U12 = mldivide(L11, H.U1 * H.V2);  %L11 U12 = A12 = H.U1 * H.V2
                L21 = mrdivide(H.U2 * H.V1, U11);  %L21 U11 = A21 = H.U2 * H.V1
                U12 = mldivide(L11, H.U1); 
                L21 = mrdivide(H.V1, U11); 
    
                NH = hrank_update(H.A22, -H.U2 * (L21 * U12), H.V2, epsilon);
                [L22, U22] = hlu(NH, 'dense', epsilon);  % lu(hadd(H.A22, L21 * U12, '-'));
                
                U12 = U12 * H.V2;
                L21 = H.U2 * L21;
    
                L = blkdiag(L11, L22);
                L(size(L11, 1)+1:end, 1:size(L11, 2)) = L21;
                L = sparse(L);
                U = blkdiag(U11, U22);
                U(1:size(U11, 1), size(U11, 2)+1:end) = U12;
                U = sparse(U);
            else
                [L, U] = lu(H.D);
            end
        else
            L = H;
            U = H;
            
            m1 = size(H.U1, 1);
            m2 = size(H.U2, 1);
            n1 = size(H.V1, 2);
    
            L.U1 = zeros(m1,1);
            L.V2 = zeros(1, m2);
            U.U2 = zeros(m2,1);
            U.V1 = zeros(1, n1);
    
            if isempty(H.D)
                [L.A11, U.A11] = hlu(H.A11, epsilon);
                
                U.U1 = htrsl(L.A11, H.U1);  %.L11 * U.U1 * U.V2 = L11 * U12 = A12 = H.U1 * H.V2
                L.V1 = htrsu(H.V1, U.A11);  % L.U2 * L.V1 * U11 = L21 * U11 = A21 = H.U2 * H.V1
                [L.A22, U.A22] = hlu(hrank_update(H.A22, -L.U2 * (L.V1 * U.U1), H.V2, epsilon), epsilon); 
    
            else
                [L.D, U.D] = lu(H.D);
                % L.D = P'*L.D;
            end
        end
    end 
    
    
    