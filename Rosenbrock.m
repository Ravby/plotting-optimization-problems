classdef Rosenbrock < Problem
    properties(Constant)
        LB = [-5 -5];  % Common bounds for Rosenbrock
        UB = [5 5];
    end
    methods
        function obj = Rosenbrock()
            obj = obj@Problem("Rosenbrock");
        end
        function out = evaluate(obj, x)
            d = length(x);
            sum = 0;
            for i = 1:d-1
                sum = sum + 100 * (x(i+1) - x(i)^2)^2 + (x(i) - 1)^2;
            end
            out = sum;
        end
    end
end


