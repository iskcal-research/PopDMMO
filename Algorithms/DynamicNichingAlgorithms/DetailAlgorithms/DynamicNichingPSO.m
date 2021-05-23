classdef DynamicNichingPSO < DynamicNichingBase
    properties(SetAccess=protected)
        pbest;
        w = 0.729;
        c1 = 2.05;
        c2 = 2.05;
    end
    
    methods
        function obj = DynamicNichingPSO(nt, rt)
            obj = obj@DynamicNichingBase(nt, rt);
        end
        
        function pop = Record(obj)
            pop = obj.pbest;
        end
        
        %% Initialize the population
        function Initialize(obj)
            obj.pop = obj.hub.GetIndis(rand_indi(obj.algRand, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper), struct("v", zeros(obj.hub.N, obj.hub.D)));
            obj.pbest = obj.pop;
        end
        
        function Evolve(obj)
%             disp(obj.hub.evaluated);
            obj.niching = obj.GetNiching();
            offspring = obj.pop;
            for i = 1:length(obj.niching)
                off = obj.GetOffspring(obj.niching{i});
                if ~isempty(off)
                    offspring(obj.niching{i}(1:length(off))) = off;
                end
            end
            obj.pop = offspring;
            obj.pbest = obj.Select(obj.pbest, obj.pop);
        end
        
         %% Generate the offspring individuals
        function offspring = GetOffspring(obj, select_idx)
            N = length(select_idx);
            pbest_decs = obj.pbest(select_idx).decs;
            pop_decs = obj.pop(select_idx).decs;
            v = obj.pop.adds("v");
            v = v(select_idx, :);
            [~, max_idx] = max(obj.pbest(select_idx).fits);
            gbest_decs = repmat(pbest_decs(max_idx, :), N, 1);
            v = obj.w .* (v + obj.c1 .* rand(obj.algRand, N, obj.hub.D) .* (pbest_decs - pop_decs) + obj.c2 .* rand(obj.algRand, N, obj.hub.D) .* (gbest_decs - pop_decs));
            off_decs = pop_decs + v;
            off_decs_new = boundary_check(off_decs, obj.hub.lower, obj.hub.upper);
            v = off_decs_new - pop_decs;            
            offspring = obj.hub.GetIndis(off_decs_new, struct("v", v));
        end
        
        %% dynamic response 
        function RespondChange(obj)
            obj.pop = obj.pbest;
            RespondChange@DynamicNichingBase(obj);
            obj.pop.SetAdds(struct("v", zeros(obj.hub.N, obj.hub.D)));
            obj.pbest = obj.pop;
        end
    end
end