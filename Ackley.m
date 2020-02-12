classdef Ackley  < Problem
    properties(Constant)
        c = 2*pi;
        b = 0.2;
        a = 20;
        LB = [-32.768 -32.768];
        UB = [32.768 32.768];
    end
    methods
        function obj = Ackley()
            obj = obj@Problem("Ackley");
        end
        function  out = evaluate(obj, x)
            d = length(x);
            sum1 = 0;
            sum2 = 0;
            for ii = 1:d
                xi = x(ii);
                sum1 = sum1 + xi^2;
                sum2 = sum2 + cos(obj.c*xi);
            end
            
            term1 = -obj.a * exp(-obj.b*sqrt(sum1/d));
            term2 = -exp(sum2/d);
            
            out = term1 + term2 + obj.a + exp(1);
        end
    end
end