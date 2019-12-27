function error = get_error(pop, bestpos, bestfit)
% used to calculate the difference between the best fit of the problem and
% of the population
    final_pop = pop{end, 3};
    final_bestpos = bestpos;
    final_bestfit = bestfit;
    error = final_bestfit - max(final_pop.fits);
end