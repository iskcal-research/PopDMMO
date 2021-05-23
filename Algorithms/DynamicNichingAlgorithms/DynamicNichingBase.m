classdef DynamicNichingBase < DynamicAlgorithm
    properties(SetAccess=protected)
        niching_type = 0;
        response_type = 0;
        memory;
        mem_size = 50;
        m = 10; % m=10
        phi = 2.0;
        niching;
    end
    methods(Access = public)
        function obj = DynamicNichingBase(nt, rt)
            obj.niching_type = nt;
            obj.response_type = rt;
        end
        
        function Evolve(obj)
%             disp([obj.hub.run, obj.hub.evaluated]);
            obj.niching = obj.GetNiching();
            next_pop = obj.pop;
            for i = 1:length(obj.niching)
                parent = obj.pop(obj.niching{i});
                offspring = obj.GetOffspring(obj.niching{i});
                next_pop(obj.niching{i}) = obj.Select(parent, offspring);
                
%                 if ~isempty(off)
%                     offspring(obj.niching{i}(1:length(off))) = off;
%                 end
            end
            obj.pop = next_pop;
%             obj.pop = obj.Select(obj.pop, offspring, obj.hub.N);
        end
        
        function niching = GetNiching(obj)
            switch(obj.niching_type)
                case 0 % no niching
                    niching = fake_niching(obj.Record().fits);
                case 1 % speciation based on the species radius
                    niching = speciation_radius(obj.Record().decs, obj.Record().fits, obj.hub.problem.gn, obj.hub.D, obj.hub.lower, obj.hub.upper);
                case 2 % speciation based on the species size
                    niching = speciation_size(obj.Record().decs, obj.Record().fits, obj.m);
                case 3 % nearest better clutering
                    niching = nearest_better_clustering(obj.Record().decs, obj.Record().fits, obj.phi);
                case 4 % nearest better clutering
                    niching = hierachical_clustering(obj.Record().decs, obj.hub.problem.gn);
                otherwise
                    error("wrong niching type in dynamic niching algorithms");
            end
        end
        
        function RespondChange(obj)
%             disp(obj.hub.evaluated);
            switch(obj.response_type)
                case 0  % no response strategy
                    obj.pop = obj.hub.GetIndis(obj.pop.decs);
                case 1  % restart
                    %  obj.pop = obj.hub.GetIndis(rand_indi(obj.algRand, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper));
                    obj.pop = restart(obj.algRand, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper);
                case 2 % reinitialize half of population
                    obj.pop = reinitialize_half_population(obj.algRand, obj.pop, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper); 
                case 3  % memory
                    % update the niching
                    obj.niching = obj.GetNiching();
                    
                    seeds_idx = [];
                    for i = 1:min(length(obj.niching), 5)
                        seeds_idx = [seeds_idx obj.niching{i}(1)];
                    end
                    [obj.pop, obj.memory] = use_memory(obj.algRand, obj.pop, obj.memory, obj.mem_size, seeds_idx, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper);
%                     obj.UseMemoryStrategy();
                otherwise
                    error("wrong response type in dynamic niching algorithms");
            end
        end
        
%         function UseMemoryStrategy(obj)
% %             disp(obj.hub.evaluated);
%             obj.pop = obj.hub.GetIndis(obj.pop.decs);
%             if ~isempty(obj.memory)
%                 obj.memory = obj.hub.GetIndis(obj.memory.decs);
%                 [~, mem_sort_idx] = sort(obj.memory.fits);
%                 obj.memory = obj.memory(mem_sort_idx);
%             end
%             
%             temp_idx = [];
%             for i = 1:min(length(obj.niching), 5)
%                 temp_idx = [temp_idx obj.niching{i}(1)];
%             end
%             
%             temp = obj.pop(temp_idx, :);
%             [~, sort_idx] = sort(obj.pop.fits, 'descend');
%             obj.pop = obj.pop(sort_idx);
%             
%             i = round(obj.hub.N / 2);
%             select_num = min(length(obj.memory), 5);
%             if ~isempty(obj.memory)
%                 obj.pop(i:i+select_num-1) = obj.memory(1:select_num);
%                 i = i + select_num-1;
%             end
%             obj.pop(i+1:end) = obj.hub.GetIndis(rand_indi(obj.algRand, obj.hub.N - i, obj.hub.D, obj.hub.lower, obj.hub.upper));
%             
%             % disp(obj.hub.evaluated);
%             for i = 1:length(temp)
%                 if length(obj.memory) < obj.mem_size
%                     obj.memory = [obj.memory; temp(i)];
%                 else
%                     dis = pdist2(temp(i).dec, obj.memory.decs);
%                     [~, nearest_idx] = min(dis);
%                     obj.memory(nearest_idx) = temp(i);
%                 end
%             end
%         end
    end
end