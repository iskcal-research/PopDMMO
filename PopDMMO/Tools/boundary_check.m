function x = boundary_check(x,lb,ub)
% function: 
%           This function is used to reflect the solutions in the decision
%           space between lower_bound and upper_bound.
% parameters:
%           x: the solutions in the decision space
%           lb: the lower boundary in the decision space
%           ub: the upper boundary in the decision space
% output:
%           x: the solutions in the decision space after reflecting the
%           boundaries

    NP = size(x, 1);
    lb = repmat(lb, NP, 1);
    ub = repmat(ub, NP, 1);

    % reflect the boundaries
    length_bound = ub - lb;
    x=(x<lb).*(lb+rem((lb-x), length_bound))+(x>=lb).*x;
    x=(x>ub).*(ub-rem((x-ub), length_bound))+(x<=ub).*x;
end

