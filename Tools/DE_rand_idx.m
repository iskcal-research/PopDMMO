function index = DE_rand_idx(algRand, N, n, s)
% function: 
%           This function is used to generate the selected DE vectors index
% parameters:
%           algRand: the random stream
%           N: the randomly generated range from 1 to N
%           n: the number of individuals
%           s: for an individual, it presents the value that randomly picked.
% output: 
%           index: it is a matrix, where each row presents randomly
%           selected indexes for DE mutation.
% example:
%           if you want to generate r1,r2 and r3 in DE/rand/1, you
%           may call DE_rand_idx(algRand, N, N, 3)
    
    index = zeros(n, s);
    for i = 1:n
        index(i, :) = randperm(algRand, N-1, s);
        index(i, index(i, :) >= i) = index(i, index(i, :) >= i) + 1;
    end
end