classdef (Abstract)DynamicProblem < Problem
    properties(SetAccess = protected)
        proRand;        % used for the dynamic random change
        env;            % the current environment number;
    end
    methods
        function obj = DynamicProblem()
            obj.proRand = RandStream.create('mt19937ar','seed', 0);
            obj.env = 0;
        end
        
        function ChangeDynamic(obj)
        end
    end
end