classdef Rastrigin < Problem
    properties(Constant)
        LB = [-5.12 -5.12];
        UB = [5.12 5.12];
    end
    methods
        function obj = Rastrigin()
            obj = obj@Problem("Rastrigin");
        end
        function out = evaluate(obj, x)
            d = length(x);
            sum = 0;
            for ii = 1:d
                xi = x(ii);
                sum = sum + (xi^2 - 10*cos(2*pi*xi));
            end
            
            out = 10*d + sum;
        end
    end
end