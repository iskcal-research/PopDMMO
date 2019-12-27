function peaks = get_peaks(pop, bestpos, bestfit)
    if isempty(bestpos)
        return;
    end
    e_d = 0.1 / 2;
    env = size(bestpos, 1);
    sum_peaks = 0;
    all_peaks =zeros(env, 3);
    for i = 1:env
        cur_gn = length(bestfit{i, 3});
        sum_peaks = sum_peaks + cur_gn;
        cur_bestfit = bestfit{i, 3}(1);
        cur_bestpos = bestpos{i, 3};
        cur_pop_decs = pop{i, 3}.decs;
        cur_pop_fits = pop{i, 3}.fits;
        for j = 1:3
            select = (cur_bestfit - cur_pop_fits) < 10^(-(j+2));
            select_pop = cur_pop_decs(select, :);
            best_dis = pdist2(cur_bestpos, select_pop);
            found = min(best_dis, [], 2) < e_d;
            if isempty(found)
                found = false(cur_gn, 1);
            end
            
            all_peaks(i,j) = sum(found);
        end
    end
    peaks = [sum(all_peaks) sum_peaks];
end