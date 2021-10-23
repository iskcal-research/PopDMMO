function decs = crossover(algRand, mutant, parent, CR)
% function: 
%           This function is used to generate the individuals by crossover
%           operator.
% parameters:
%           algRand: the random stream
%           mutant: the mutant in the decision space
%           parent: the parent in the decision space
%           CR: the crossover rate
% output: 
%           decs: the individuals after crossover operating
    [N, D] = size(parent);
    cross = rand(algRand, N, D) < CR;
    
    parent_only_idx = find(sum(cross, 2) == 0);
    for i = parent_only_idx
        cross(i, randi(algRand, D)) = true;
    end
    decs = cross .* mutant + (1-cross) .* parent;
end