function [pop, mem] = use_memory(algRand, pop, mem, mem_size, seeds_idx, N, D, lower, upper)
    pop = Hub.GetIns().GetIndis(pop.decs);
    if ~isempty(mem)
        mem = Hub.GetIns().GetIndis(mem.decs);
        [~, mem_sort_idx] = sort(mem.fits);
        mem = mem(mem_sort_idx);
    end

%     temp_idx = [];
%     for i = 1:min(length(obj.niching), 5)
%         temp_idx = [temp_idx obj.niching{i}(1)];
%     end

    temp = pop(seeds_idx);
    [~, sort_idx] = sort(pop.fits, 'descend');
    pop = pop(sort_idx);

    i = round(N / 2);
    select_num = min([length(seeds_idx), 5, length(mem)]);
    if ~isempty(mem)
        pop(i+1:i+select_num) = mem(1:select_num);
        i = i + select_num;
    end
    pop(i+1:end) = Hub.GetIns().GetIndis(rand_indi(algRand, N - i, D, lower, upper));

    % disp(obj.global_.evaluated);
    for i = 1:length(temp)
        if length(mem) < mem_size
            mem = [mem; temp(i)];
        else
            dis = pdist2(temp(i).dec, mem.decs);
            [~, nearest_idx] = min(dis);
            if temp(i).fit > mem(nearest_idx).fit
                mem(nearest_idx) = temp(i);
            end
        end
    end
end