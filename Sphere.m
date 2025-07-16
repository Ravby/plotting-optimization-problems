classdef Sphere  < Problem
    properties(Constant)
        LB = [-100 -100];
        UB = [100 100];
    end
    methods
        function obj = Sphere()
            obj = obj@Problem("Sphere");
        end
        function out = evaluate(obj, x)
            d = length(x);
            sum = 0;
            for ii = 1:d
                xi = x(ii);
                sum = sum + xi^2;
            end
            
            out = sum;
        end
    end
end