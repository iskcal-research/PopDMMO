function idx = tourment_select(algRand, fit, n, k)
% function: 
%           This function implements the selection operator in CEP
% parameters:
%           algRand: the random stream
%           fit: the finess values of the population
%           n: the number of the selected individuals
%           k: the number of the individuals in a single competition
% output:
%           idx: the index of the selected individuals
    s = size(fit, 1);
    k = min(s, k);
    
    select_fit = zeros(s, k);
    select = zeros(s, k);
    for i = 1:s
       select(i, :) = randperm(algRand, s, k);
       select_fit(i, :) = fit(select(i, :));
    end
    win = sum(repmat(fit, 1, k) >= select_fit, 2);   
    [~, idx] = sortrows([win fit], 'descend');
    
    idx = idx(1:n);
end