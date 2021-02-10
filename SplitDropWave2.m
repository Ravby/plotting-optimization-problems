classdef SplitDropWave2  < Problem
    properties(Constant)
        LB = [-5 -5];
        UB = [5 5];
    end
    methods
        function obj = SplitDropWave2()
            obj = obj@Problem("SplitDropWave2");
        end
        function out = evaluate(obj, x)
            x1 = x(1);
            x2 = x(2);
            
            out = cos(x1^2 + x2^2) + 2 * exp(-10 * x2^2) + 0.5 * x2;
        end
    end
end
