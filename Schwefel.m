classdef Schwefel  < Problem
    properties(Constant)
        LB = [-500 -500];
        UB = [500 500];
    end
    methods
        function obj = Schwefel()
            obj = obj@Problem("Schwefel");
        end
        function out = evaluate(obj, x)
            d = length(x);
            sum = 0;
            for ii = 1:d
                xi = x(ii);
                sum = sum + xi*sin(sqrt(abs(xi)));
            end
            
            out = 418.9829*d - sum;
        end
    end
end