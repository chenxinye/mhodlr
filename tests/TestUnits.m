classdef TestUnits < matlab.unittest.TestCase
    methods (Test)

        function testBuild(testCase)

            rng(0)
            A = rand(100, 100);
            u1 = precision('d');
            u2 = precision('s');
            u3 = precision('h');
            u4 = precision('b');
            u5 = precision('q52');
            
            u_chain = prec_chain(u1, u2, u3, u4, u5);
            
            epsilon = 1e-4;
            depth = 3;
            
            hA = hodlr(A, depth, 10, 'svd', epsilon); 
            rA = recover(hA);
            err = norm(rA - A, 'fro');
            load('error/rerr1.mat');
            testCase.assertNotEqual(err, err1);
            
            aphA = amphodlr(u_chain, A, depth, 2, 'svd', epsilon); 
            aprA = recover(aphA);
            err = norm(aprA - A, 'fro');
            load('error/rerr2.mat');
            testCase.assertNotEqual(err, err2);
            
            mphA = mphodlr(u_chain, A, depth, 2, 'svd', epsilon); 
            mprA = recover(mphA);
            err = norm(mprA - A, 'fro');
            load('error/rerr3.mat');
            testCase.assertNotEqual(err, err3);
        end
    end

end

