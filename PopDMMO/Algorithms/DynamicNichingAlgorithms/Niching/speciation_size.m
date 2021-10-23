function niching = speciation_size(decs, fit, m)
    [~, sort_idx] = sort(fit, 'descend');
    selected = false(size(decs, 1), 1);
    niching = cell(0);
    dis = pdist2(decs, decs);
    
    while ~all(selected)
        idx = sort_idx(selected == false);
%         niching{end+1} = idx(1);
%         sort_idx(idx(1)) = true;
        cur_dis = dis(idx(1), idx);
        if length(cur_dis) <= m
            cur_select = idx;
        else
            [~, mink_cur_idx] = sort(cur_dis);
            cur_select = idx(mink_cur_idx(1:m));
        end
        [~, pos] = ismember(cur_select, sort_idx);
        selected(pos) = true;
        niching{end+1} = cur_select;
    end
end