function idx = roulette_select(algRand, fits, n)
% function: 
%           This function implements the roulette selection operator
% parameters:
%           algRand: the random stream
%           fits: the finess value of the population
%           n: the number of the selected individuals
% output:
%           idx: the index of the selected individuals
    N = size(fits, 1);
    norm_fit = 1 ./ (1+max(fits)-fits);
    phi_s = norm_fit ./ sum(norm_fit);
    
    select_pro = cumsum(phi_s);
    select_pro = repmat(select_pro', n, 1);
    
    r = rand(algRand, n, 1);
    r = repmat(r, 1, N);
    
    idx = sum(select_pro < r, 2) + 1;
end