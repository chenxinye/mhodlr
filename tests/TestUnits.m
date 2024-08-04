classdef TestUnits < matlab.unittest.TestCase
    methods (Test)

        function testBuild(testCase)
            A = rand(100, 100);
            u1 = precision('d');
            u2 = precision('s');
            u3 = precision('h');
            u4 = precision('b');
            u5 = precision('q52');
            
            u_chain = prec_chain(u1, u2, u3, u4, u5);

            epsilon = 1e-4;
            depth = 3;
            aphA = amphodlr(u_chain, A, depth, 10, 'svd', epsilon); 
            aprA = recover(aphA);
        end

        
    end

end

