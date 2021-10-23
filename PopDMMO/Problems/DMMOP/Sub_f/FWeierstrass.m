function f = FWeierstrass(x)
    [NP, D] = size(x);
    a = 0.5;
    b = 3;
    k_max = 20;
    c1(1:k_max+1) = a .^ (0:k_max);
    c2(1:k_max+1) = 2 .* pi .* b .^ (0:k_max);
    t = zeros(NP, D);
    for i = 1:D
       t(:, i) = sum(c1 .* cos(c2 .* (x(:, i) + 0.5)), 2); 
    end
    f = sum(t, 2) - D .* sum(c1 .* cos(c2 .* 0.5), 2);
end