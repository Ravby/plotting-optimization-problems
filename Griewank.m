classdef Griewank  < Problem
    properties(Constant)
        LB = [-600 -600];
        UB = [600 600];
    end
    methods
        function obj = Griewank()
            obj = obj@Problem("Griewank");
        end
        function out = evaluate(obj, x)
            n = size(x, 2);
    
            sumcomp = 0;
            prodcomp = 1;

            for i = 1:n
                sumcomp = sumcomp + (x(:, i) .^ 2);
                prodcomp = prodcomp .* (cos(x(:, i) / sqrt(i)));
            end

            out = (sumcomp / 4000) - prodcomp + 1; 
        end
    end
end
