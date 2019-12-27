classdef (Abstract)DynamicAlgorithm < Algorithm
    methods
        function changed = CheckChange(obj)
            changed = rem(obj.hub.evaluated, obj.hub.freq) == 0;
        end
        
        function RespondChange(obj)
            obj.pop = obj.hub.GetIndis(obj.pop.decs);
        end
    end
end