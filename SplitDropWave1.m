classdef SplitDropWave1  < Problem
    properties(Constant)
        LB = [-3 -3];
        UB = [3 3];
    end
    methods
        function obj = SplitDropWave1()
            obj = obj@Problem("SplitDropWave1");
        end
        function out = evaluate(obj, x)
            
            x1 = x(1);
            x2 = x(2);
            
            out = cos(x1^2 + x2^2) + 2 * exp(-10 * x2^2);
        end
    end
end
