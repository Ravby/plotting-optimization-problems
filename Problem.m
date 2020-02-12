classdef (Abstract) Problem

    properties
        problem_name = "";
    end
    
    methods
        function obj = Problem(name)
            obj.problem_name = name;
        end
    end
    methods (Abstract)
        out = evaluate(obj, x)
    end
end

