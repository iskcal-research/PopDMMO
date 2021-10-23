classdef (Abstract)Problem < handle
    properties(SetAccess=protected)
        hub;        % the hub object
    end    
    methods        
        function obj = Problem()
            obj.hub = Hub.GetIns();
            
            obj.hub.D = 2;
            obj.hub.encoding = "real";

            obj.hub.lower = 0;
            obj.hub.upper = 1;
        end
        
        function Initialize(obj)
        end
        
        function fits = GetFits(obj, decs)
            fits = zeros(size(decs, 1), 1);
        end
        
        function cons = GetCons(obj, decs)
            cons = zeros(size(decs, 1), 1);
        end
        
        function bestfit = GetBestFit(obj)
            bestfit = 0;
        end
        
        function bestpos = GetBestPos(obj)
            bestpos = zeros(1, obj.hub.D);
        end
    end
end