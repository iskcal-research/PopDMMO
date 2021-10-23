function f = FRastrigin(x)
    f = sum(x.^2 - 10 .* cos(2 .* pi .* x) + 10 , 2);
end