classdef DynamicNichingDE < DynamicNichingBase
    methods
        function obj = DynamicNichingDE(nt, rt)
            obj = obj@DynamicNichingBase(nt, rt);
        end
        
        function offspring = GetOffspring(obj, select_idx)
            if length(select_idx) < 4
                offspring = obj.pop(select_idx);
            else
                N = length(select_idx);
                idx = DE_rand_idx(obj.algRand, N, N, 3);  % DE/rand/1
                pop_decs = obj.pop(select_idx).decs;

                off_decs = pop_decs(idx(:, 1), :) + 0.5 * (pop_decs(idx(:, 2), :) - pop_decs(idx(:, 3), :));
                off_decs = boundary_check(off_decs, obj.hub.lower, obj.hub.upper);

                off_decs = crossover(obj.algRand, off_decs, pop_decs, 0.9);
                offspring = obj.hub.GetIndis(off_decs);
            end
        end
    end
end