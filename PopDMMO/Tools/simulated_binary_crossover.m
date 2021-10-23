function decs = simulated_binary_crossover(algRand, pop, fit, lower, upper, eta)
% function: 
%           This function implements the binary crossover method to
%           exchange two individuals
% parameters:
%           algRand: the random stream
%           pop: the decision value of the population
%           fit: the finess value of the population
%           lower: the lower boundary of the decision space
%           upper: the upper boundary of the decison space
%           eta: the parameter used in the simulated binary crossover operator
% output:
%           decs: the decision after the simulated binary crossover operator

    [N, D] = size(pop);
    idx1 = roulette_select(algRand, fit, ceil(N/2));
    idx2 = roulette_select(algRand, fit, ceil(N/2));
    parent1 = pop(idx1, :);
    parent2 = pop(idx2, :);
    
%     parent1 = pop(1:2:end, :);
%     parent2 = pop(2:2:end, :);
%     if size(parent1, 1) > size(parent2, 1)
%         parent2 = [parent2; pop(1, :)];
%     end
    
    mu = rand(algRand, ceil(N/2), D);
    r = (mu <= 0.5) .* (2 .* mu).^(1./(eta+1)) + (mu > 0.5) .* (1./(2.*(1-mu))).^(1./(eta+1));
    decs1 = 0.5 .* ((1+r).*parent1 + (1-r).*parent2);
    decs2 = 0.5 .* ((1-r).*parent1 + (1+r).*parent2);
    
    decs = [decs1; decs2];
%    decs = zeros(N, D);
%     decs(1:2:end, :) = decs1;
%     decs(2:2:end, :) = decs2(1:floor(N/2), :);
    decs = boundary_check(decs, lower, upper);
    decs = decs(1:N, :);
end