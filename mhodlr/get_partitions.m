function [m, n] = get_partitions(hA)
       
    if ~isempty(hA.D)
        [m, n] = hsize(hA);
    else
        [m1, n1] = get_partitions(hA.A11);
        [m2, n2] = get_partitions(hA.A22);
        
        m = max(length(m1), length(m2));
        
        m1 = [ m1, ones(1, m - length(m1)) * m1(end) ];
        m2 = [ m2, ones(1, m - length(m2)) * m2(end) ];
        n1 = [ n1, ones(1, m - length(n1)) * n1(end) ];
        n2 = [ n2, ones(1, m - length(n2)) * n2(end) ];
        
        m = [ m1, m2 + m1(end) ];
        n = [ n1, n2 + n1(end) ];
    end
    
end