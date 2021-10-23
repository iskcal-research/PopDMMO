function pop = rand_indi(algRand, N, D, lb, ub)
% function: 
%           This function is used to randomly generate a 
%           population in the decision space
% parameters:
%           algRand: the random stream
%           N: the population size
%           D: the dimension of the problem
%           lb: the lower boundary of the decision space
%           ub: the upper boundary of the decison space
% output:
%           pop: the population matrix in the decision space
   
    pop = (ub-lb) .* rand(algRand, N, D) + lb;
end