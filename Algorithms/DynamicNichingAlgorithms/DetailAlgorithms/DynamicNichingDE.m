classdef DynamicNichingDE < DynamicNichingBase
    methods
        function obj = DynamicNichingDE(nt, rt)
            obj = obj@DynamicNichingBase(nt, rt);
        end
        
        function offspring = GetOffspring(obj, select_idx)
            N = length(select_idx);
            num = N;
            if N < 4
                num = 4;
                seed = select_idx(1);
                other_idx = setdiff(1:obj.hub.N, select_idx);
                seed_dis = pdist2(obj.pop(seed).decs, obj.pop(other_idx).decs);
                [~, near_idx_idx] = mink(seed_dis, num-N);
                select_idx = [select_idx; other_idx(near_idx_idx)'];
            end

            idx = DE_rand_idx(obj.algRand, num, N, 3);  % DE/rand/1
            pop_decs = obj.pop(select_idx).decs;

            off_decs = pop_decs(idx(:, 1), :) + 0.5 * (pop_decs(idx(:, 2), :) - pop_decs(idx(:, 3), :));
            off_decs = boundary_check(off_decs, obj.hub.lower, obj.hub.upper);

            off_decs = crossover(obj.algRand, off_decs, pop_decs(1:N, :), 0.9);
            offspring = obj.hub.GetIndis(off_decs);
        end
    end
end