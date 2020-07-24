classdef HolderTable  < Problem
    properties(Constant)
        LB = [-10 -10];
        UB = [10 10];
    end
    methods
        function obj = HolderTable()
            obj = obj@Problem("HolderTable");
        end
        function out = evaluate(obj, x)
            n = size(x, 2);
            X = x(:, 1);
            Y = x(:, 2);

            expcomponent = exp( abs(1 - (sqrt(X .^2 + Y .^ 2) / pi)) );

            out = -abs(sin(X) .* cos(Y) .* expcomponent);
        end
    end
end
