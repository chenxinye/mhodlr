classdef TestUnits < matlab.unittest.TestCase
    % TestExamples contains a set of 4 simple tests:
    %     1) an equality test for a non-leap year date
    %     2) an equality test for a leap year date
    %     3) a negative test for an invalid date format input
    %     4) a negative test for a correct date format but an invalid date
    %     5) an equality test for a non-leap year date using the alternate
    %        dateFormat (COMMENTED OUT)
    %
    % Notes:
    %     A) A negative test verifies that the code errors/fails in an
    %        expected way (e.g., the code gives the right error for a
    %        specific bad input)
    %     B) The 5th test is included for completeness, but is commented 
    %        out to illustrate missing code coverage in continous
    %        integration (CI) systems

    % Copyright 2022 The MathWorks, Inc.

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

