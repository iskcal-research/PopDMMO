classdef DynamicNichingGA < DynamicNichingBase
    properties(SetAccess=protected)
        eta = 1;
        mut = 0.01;
    end
    methods(Access=public)
        function obj = DynamicNichingGA(nt, rt)
            obj = obj@DynamicNichingBase(nt, rt);
        end
        
        function offspring = GetOffspring(obj, select_idx)
            pop_decs = obj.pop(select_idx).decs;
            pop_fits = obj.pop(select_idx).fits;
            off_decs = simulated_binary_crossover(obj.algRand, pop_decs, pop_fits, obj.hub.lower, obj.hub.upper, obj.eta);
            mut_per = rand(obj.algRand, size(off_decs));
            mut_pos = mut_per < obj.mut;
            off_decs(mut_pos) = rand(obj.algRand, sum(mut_pos(:)), 1) .* (obj.hub.upper(1) - obj.hub.lower(1)) + obj.hub.lower(1);
            offspring = obj.hub.GetIndis(off_decs);
        end
        
        function next_pop = Select(obj, parent, offspring)
            N = size(parent, 1);
            pop = [parent; offspring];
            [~, sort_idx] = sort(pop.fits, 'descend');
            next_pop = pop(sort_idx(1:N));
            
%             if length(offspring) < length(parent)
%                 offspring = [offspring; parent(length(offspring)+1:end)];
%             end
%             
%             next_pop = parent;
%             for i= [1:2:length(offspring)]
%                 end_idx = min(i+1, length(parent));
%                 pop = [parent(i:end_idx);offspring(i:end_idx)];
%                 [~, sort_idx] = sort(pop.fits, 'descend');
%                 next_pop(i:end_idx) = pop(sort_idx(end_idx-i+1));
%             end
        end
    end
end