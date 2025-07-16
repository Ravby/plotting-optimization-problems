classdef Schaffer < Problem
    properties(Constant)
        LB = [-100 -100];
        UB = [100 100];
    end
    methods
        function obj = Schaffer()
            obj = obj@Problem("Schaffer");
        end
        function out = evaluate(obj, x)
            if length(x) ~= 2
                error('Schaffer function is defined for two variables (x, y).');
            end
            x1 = x(1);
            x2 = x(2);
            numerator = sin(x1^2 - x2^2)^2 - 0.5;
            denominator = (1 + 0.001 * (x1^2 + x2^2))^2;
            out = 0.5 + numerator / denominator;
        end
    end
end


