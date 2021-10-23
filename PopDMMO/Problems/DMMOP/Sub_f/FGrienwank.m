function f = FGrienwank(x)
    D = size(x, 2);
    f1 = sum(x.^2/4000 , 2);
    f2 = prod(cos(x ./ sqrt(1:D)), 2);
    f = f1 - f2 + 1;
end