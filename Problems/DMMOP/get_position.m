function pos = get_position(N, D, lb, ub, dpeak)
    seed = RandStream.create('mt19937ar','seed',0);
    pos = zeros(N, D);
    unsatisfied_index = 1:N;
    while ~isempty(unsatisfied_index)
        pos(unsatisfied_index, :) = rand(seed, length(unsatisfied_index), D) .* (ub - lb) + lb;
        pos_dis = pdist2(pos, pos);
        unsatisfied_index = find(all(pos_dis > dpeak, 2));
    end
end