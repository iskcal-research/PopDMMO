function pop = reinitialize_half_population(algRand, pop, N, D, lower, upper)
    pop = Hub.GetIns().GetIndis(pop.decs);
    
    [~, sort_idx] = sort(pop.fits, 'descend');
    
    pop = pop(sort_idx);
    
    pop(ceil(N/2)+1:end) = Hub.GetIns().GetIndis(rand_indi(algRand, N-ceil(N/2), D, lower, upper));
end