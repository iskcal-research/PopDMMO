classdef DynamicNichingCSA < DynamicNichingBase
    properties(SetAccess=protected)
        d = 0.1;
        copy = 2;
        % nc = 0.1;
    end
    
    methods(Access=public)   
        function obj = DynamicNichingCSA(nt, rt)
            obj = obj@DynamicNichingBase(nt, rt);
        end
        
        function Evolve(obj)
            obj.niching = obj.GetNiching();
            next_pop = obj.pop;
            for i = 1:length(obj.niching)
                parent = obj.pop(obj.niching{i});
                offspring = obj.GetOffspring(obj.niching{i});
                next_pop(obj.niching{i}) = obj.Select(parent, offspring);
               %offspring = [offspring; off];
            end
            % obj.pop = obj.Select(obj.pop, offspring, obj.hub.N);
             %% several worst individuals would be replaced by the randomly generated individuals
            [~, sort_idx] = sort(next_pop.fits, 'descend'); 
            start = round(obj.hub.N * (1-obj.d));
            receptor = obj.hub.GetIndis(rand_indi(obj.algRand, obj.hub.N - start, obj.hub.D, obj.hub.lower, obj.hub.upper));
            if ~isempty(receptor)
%                 obj.pop(start+1:start+size(receptor, 1)) = receptor;
                next_pop(sort_idx(obj.hub.N - length(receptor) + 1:end)) = receptor;
            end
            obj.pop = next_pop;
        end
        
         %% Generate the offspring individuals with the selected parents
        function offspring = GetOffspring(obj, select_idx)
            n = length(select_idx);
            if n == 1
                %offspring = [];
                offspring = obj.pop(select_idx);
                return;
            end
            pop = obj.pop(select_idx);
            pop_decs = pop.decs;
            pop_fits = pop.fits;
            norm_pop_fits = (mapminmax(pop_fits', 0, 1))';
            alpha = exp(-4 * norm_pop_fits);
            
            offspring = pop;
            for i = 1:n
%                 disp(obj.hub.evaluated);
                cur_n = ceil(0.1*n/i);
                cur_decs = repmat(pop_decs(i, :), cur_n, 1);
                select_idx = randi(obj.algRand, n-1, cur_n, 1);
                select_idx(select_idx>=i) = select_idx(select_idx>=i) + 1;
                beta = rand(obj.algRand, size(cur_decs));
                temp_decs = (1-beta) .* cur_decs + beta .* pop_decs(select_idx, :);
                change_pos = rand(obj.algRand, size(cur_decs)) < alpha(i);
                final_decs = change_pos .* temp_decs + (1-change_pos) .* cur_decs;
                 
                off = obj.hub.GetIndis(final_decs);
                if ~isempty(off)
                    [~, max_idx] = max(off.fits);
                    offspring(i) = off(max_idx);
                end
            end
        end
        
%         function next_pop = Select(obj, parent, offspring)
%             N = size(parent, 1);
%             pop = [parent; offspring];
%             [~, sort_idx] = sort(pop.fits, 'descend');
%             next_pop = pop(sort_idx(1:N));
%         end
        
%         function next_pop = Select(obj, pop, offspring, N)
%             pop = [pop; offspring];
%             pop_fits = pop.fits;
%             [~, sort_idx] = sort(pop_fits, 'descend');
%             next_pop = pop(sort_idx(1:N));
%         end
    end
end