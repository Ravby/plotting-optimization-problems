classdef Easom  < Problem
    properties(Constant)
        LB = [-100 -100];
        UB = [100 100];
    end
    methods
        function obj = Easom()
            obj = obj@Problem("Easom");
        end
        function out = evaluate(obj, x)
            x1 = x(1);
            x2 = x(2);
            
            fact1 = -cos(x1)*cos(x2);
            fact2 = exp(-(x1-pi)^2-(x2-pi)^2);
            
            out = fact1*fact2;
        end
    end
end