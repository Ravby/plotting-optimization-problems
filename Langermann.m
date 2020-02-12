classdef Langermann  < Problem
    properties(Constant)
        LB = [0 0];
        UB = [10 10];
        m = 5;
        c = [1, 2, 5, 2, 3];
        A = [3, 5; 5, 2; 2, 1; 1, 4; 7, 9];
    end
    methods
        function obj = Langermann()
            obj = obj@Problem("Langermann");
        end
        function out = evaluate(obj, x)
            d = length(x);
            outer = 0;
            for ii = 1:obj.m
                inner = 0;
                for jj = 1:d
                    xj = x(jj);
                    Aij = obj.A(ii,jj);
                    inner = inner + (xj-Aij)^2;
                end
                new = obj.c(ii) * exp(-inner/pi) * cos(pi*inner);
                outer = outer + new;
            end
            
            out = outer;
        end
    end
end