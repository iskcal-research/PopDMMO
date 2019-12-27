function niching = speciation_radius(decs, fit, N_o, D, lower, upper)
    N = size(decs, 1);
    [~, sort_idx] = sort(fit, 'descend');
    
    seeds = [];
    species = zeros(N, 1);
    niching = cell(0);
    
    r_s = (pdist2(upper, lower) / 2 ) / (N_o ^ (1/D));
    
    for i = 1:N
        cur_idx = sort_idx(i);
        if isempty(seeds)
            seeds = [seeds; cur_idx];
            niching{end+1} = cur_idx;
            species(cur_idx) = length(niching);
        else
            dis = pdist2(decs(cur_idx, :), decs(seeds, :));
            if all(dis > r_s) 
                seeds = [seeds; cur_idx];
                niching{end+1} = cur_idx;
                species(cur_idx) = length(niching);
            else
                [~, min_idx] = min(dis);
                niching{species(seeds(min_idx))} = [niching{species(seeds(min_idx))}; cur_idx];
                species(cur_idx) = species(seeds(min_idx));
            end
        end
    end
end