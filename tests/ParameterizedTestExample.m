classdef ParameterizedTestExample < matlab.unittest.TestCase
    % Creates 12 test points, one test point for the 15th day of every month of 2021
    
	% Copyright 2022 The MathWorks, Inc.
	
    properties (TestParameter)
        monthNum = num2cell(1:12);
        dayNum = {15};
        yearNum = {2021};
    end
    
    methods (Test)
        function testDayofyear(testCase,monthNum,dayNum,yearNum)
            
        end
    end
    
end

