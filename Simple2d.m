classdef Simple2d  < Problem
    properties(Constant)
        LB = [-5 -5];
        UB = [5 5];
    end
    methods
        function obj = Simple2d()
            obj = obj@Problem("Simple2d");
        end
        function out = evaluate(obj, x)
            out = cos(x(1)) +sin(x(2) / 4);
        end
    end
end

