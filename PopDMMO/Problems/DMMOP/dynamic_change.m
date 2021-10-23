function u = dynamic_change(proRand, u, change_type, u_min, u_max, u_severity, change_count)
    % parameters settings
    p = 12;
    noisy_servity = 0.8;
    A = 3.67;
    alpha = 0.04;
    alpha_max = 0.1;
    
    u_range = u_max - u_min;
    switch(change_type)
        case 1 % small step
            r = 2 * rand(proRand, size(u, 1), 1) - 1;
            u = max(min(u+alpha .* u_range .* r .* u_severity, u_max), u_min);
        case 2 % large step
            r = 2 * rand(proRand, size(u, 1), 1) - 1;
            %u = max(min(u+u_range*(alpha*sign(r)+(alpha_max-alpha)*r)*u_severity,u_max),u_min);
            u = max(min(u+u_range.*(alpha.*sign(r)+(alpha_max-alpha).*r).*u_severity ,u_max), u_min);
        case 3 % random
            %u = max(min(u+normrnd(0,1,size(u,1),1)*u_severity,u_max),u_min);
            u = max(min(u+randn(proRand, size(u, 1), 1).*u_severity, u_max), u_min);
        case 4 % chaotic
            u = (u - u_min)./u_range;
            u = A.*u.*(1-u);
            u = max(min(u_min + u.*u_range, u_max), u_min);
%             u = u_min + u.*u_range;
        case 5 % recurrent
            load phi;
            phi = phi(1:size(u, 1));
            u = u_min + u_range .* (sin(2*pi*change_count/p+phi)+1)/2;
        case 6 % recurrent with noisy
            load phi;
            phi = phi(1:size(u, 1));
            u = u_min + u_range .* (sin(2*pi*change_count/p+phi)+1)/2 + randn(proRand, size(u, 1), 1)*noisy_servity;
            u = max(min(u, u_max), u_min);
        otherwise
            error('change_type has the wrong number.');
    end
end