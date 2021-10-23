classdef DMMOP < DynamicProblem
    properties(SetAccess=private)
        fun_num = 1;        % function number, F_1~F_8
        change_type = 1;    % dynamic change mode, C1~C8 
        o;                  % the best position in the current environment
        x;                  % the best position of all peaks
        x_initial;          % the initial best position of all peaks
        x_angel = 0;        % the rotation angel of x
%         o_initial;
        of;                 % the best fitness
        gn;                 % the number of global peaks in the current environment
        gn_max = 4;         % the maximal number of global peaks
        dpeaks = 0.1;       % minimal distance between each two peaks
        selected_idx;       % the selected index of the global peaks
        is_add = false;     % used in C7, deciding the number of global peaks is increasing or decreasing
       
        % the following four parameters are used in DF1
        ln = 0;             % the number of local peaks
        h;                  % the heights of peaks
        w;                  % the widths of peaks
        h_s = 7.0;          % the severity of heights
        w_s = 1.0;          % the severity of widths
        h_bound = [30, 70]; % the range of heights
        w_bound = [1, 12];  % the range of widths
        
        % the following one parameter is used in composition functions
        M;                  % the rotation matrix
        M_initial;          % the initial rotation matrix          
        M_angel;            % the rotation angel of M
        sub_fun;            % the list of subfunctions
        lambda;             % parameter lambda
        sigma;              % parameter sigma
        rs;                  % random dimension in C5 and C6
        
        % history best position and fitness;
        his_o = {};         % the history best position of the optima
        his_of = {};        % the history best fitness of the optima
    end
    methods       
        function obj = DMMOP(fun_num, change_type, D)
            obj.fun_num = fun_num;
            obj.change_type = change_type;
            obj.hub.D = D;
            obj.hub.lower = -5 .* ones(1, D);
            obj.hub.upper = 5 .* ones(1, D);
            obj.hub.freq = 5000 * D;
            obj.hub.evaluation = obj.maxEnv * obj.hub.freq;
            obj.rs = struct();
        end
        
        function Initialize(obj, evaluated)
            %% set the parameters in DMMOP
            if obj.fun_num <= 4 % DMMOP is based on DF1
                
                %% set parameters of DF1
                h_global = 75;
                switch(obj.fun_num)
                    case 1 % absolutely random 
                        obj.ln = 4;
                        obj.h = [ones(1, obj.gn_max) .* h_global   62.5889   66.2317   68.0795   66.5350]';
                        obj.w = [7.9560    2.0729    4.0635    7.0157   11.5326   11.6138    2.7337   11.6765]';
                        file = fullfile(fileparts(mfilename('fullpath')),'./Data/F1X');
                        t = load(file);
                        obj.x = t.data;
                        obj.x = obj.x(1:obj.gn_max+obj.ln, 1:obj.hub.D);
                        
%                         obj.h = [ones(obj.gn, 1) .* h_global; rand(obj.proRand, obj.ln, 1) .* (obj.h_bound(2) - obj.h_bound(1)) + obj.h_bound(1)];
%                         obj.w = rand(obj.proRand, obj.gn+obj.ln, 1) .* (obj.w_bound(2) - obj.w_bound(1)) + obj.w_bound(1);
%                         obj.x = rand(obj.proRand, obj.gn+obj.ln, obj.hub.D) .* (obj.hub.upper - obj.hub.lower) + obj.hub.lower;
                    case 2 % simulate Shubert function 
                        if rem(obj.gn_max, 2) == 1
                            obj.gn_max = obj.gn_max + 1;
                        end
                        obj.ln = 0;
                        
                        obj.h = ones(obj.gn_max, 1) .* h_global;
                        obj.w = ones(obj.gn_max+obj.ln, 1) .* 12;
                        obj.x = repmat([-3;-2;2;3], 1, obj.hub.D);
                    case 3 % simulate Vincent function
                        obj.ln = 0;
                        obj.h = ones(obj.gn_max, 1) .* h_global;
                        obj.w = ones(obj.gn_max+obj.ln, 1) .* 5;
                        obj.x = repmat([-3;-2;0;4]+0.5, 1, obj.hub.D);
                    case 4 % simulate modified Rastrigin function
                        obj.ln = 0;
                        obj.h = ones(obj.gn_max, 1) .* h_global;
                        obj.w = ones(obj.gn_max, 1) .* 5;
                        obj.x = repmat([-3;-1;1;3], 1, obj.hub.D);
                    otherwise
                        error('fun_num in DMMOP is illegal.');
                end
                % check whether the distance satisfies the minimal distance
                obj.x = obj.SetMinDistance(obj.x);
                obj.x_initial = obj.x; % save the initial states
                % if the change mode is C5 or C6, change the environment first
                if obj.change_type == 5 || obj.change_type == 6
                    obj.h = [ones(obj.gn_max, 1) .* h_global; dynamic_change(obj.proRand, obj.h(obj.gn_max+1:end), obj.change_type, obj.h_bound(1), obj.h_bound(2), obj.h_s, obj.env)];
                    obj.w = dynamic_change(obj.proRand, obj.w, obj.change_type, obj.w_bound(1), obj.w_bound(2), obj.w_s, obj.env);                    
                    [obj.x, obj.x_angel] = obj.ChangeMatrix(obj.x, obj.x_initial, obj.hub.lower, obj.hub.upper, obj.x_angel, 'x');
                end
                obj.x = obj.SetMinDistance(obj.x);
%                 obj.of = obj.h(1:obj.gn_max);
%                 obj.o = obj.x(1:obj.gn_max, :);
            else % DMMOP is based on the composition functions
                obj.ln = 0;
                switch(obj.fun_num)
                    case 5  % CF1
                        obj.sub_fun = {@FGrienwank, @FGrienwank, @FWeierstrass, @FWeierstrass, @FSphere, @FSphere};
                        obj.gn_max = length(obj.sub_fun);
                        obj.sigma = ones(1, obj.gn_max);
                        obj.lambda = [1, 1, 8, 8, 1/5, 1/5]';
                        obj.M = cell(obj.gn_max, 1);
                        obj.M(1:end) = {diag(ones(1, obj.hub.D))};
                    case 6 % CF2
                        obj.sub_fun = {@FRastrigin, @FRastrigin, @FWeierstrass, @FWeierstrass, @FGrienwank, @FGrienwank, @FSphere, @FSphere};
                        obj.gn_max = length(obj.sub_fun);
                        obj.sigma = ones(1, obj.gn_max);
                        obj.lambda = [1, 1, 10, 10, 1/10, 1/10, 1/7, 1/7]';
                        obj.M = cell(obj.gn_max, 1);
                        obj.M(1:end) = {diag(ones(1, obj.hub.D))};
                    case 7 %CF3
                        obj.sub_fun = {@FEF8F2, @FEF8F2, @FWeierstrass, @FWeierstrass, @FGrienwank, @FGrienwank};
                        obj.gn_max = length(obj.sub_fun);
                        obj.sigma = [1,1,2,2,2,2];
                        obj.lambda = [1/4, 1/10, 2, 1, 2, 5]';
                        if ~ismember(obj.hub.D, [2, 3, 5, 10, 20])
                            error('the dimension is illegal in function 7.');
                        end
                        file = fullfile(fileparts(mfilename('fullpath')), sprintf('./Data/CF3_M_D%d', obj.hub.D));
                        t = load(file);
                        obj.M = t.M;
                    case 8 %CF4
                        obj.sub_fun = {@FRastrigin, @FRastrigin, @FEF8F2, @FEF8F2, @FWeierstrass, @FWeierstrass, @FGrienwank, @FGrienwank};
                        obj.gn_max = length(obj.sub_fun);
                        obj.sigma = [1,1,1,1,1,2,2,2];
                        obj.lambda = [4, 1, 4, 1, 1/10, 1/5, 1/10, 1/40]';
                        if ~ismember(obj.hub.D, [2, 3, 5, 10, 20])
                            error('the dimension is illegal in function 8.');
                        end
                        file = fullfile(fileparts(mfilename('fullpath')), sprintf('./Data/CF4_M_D%d', obj.hub.D));
                        t = load(file);
                        obj.M = t.M;
                    otherwise
                        error('fun_num in DMMOP is illegal.');
                end
                file = fullfile(fileparts(mfilename('fullpath')),'./Data/optima');
                t = load(file);
                
%                 obj.o = t.o(1:obj.gn, 1:obj.hub.D);
                obj.x = t.o(end:-1:end-obj.gn_max+1, 1:obj.hub.D);
                obj.x = obj.SetMinDistance(obj.x);
                obj.x_initial = obj.x;
                
                obj.M_initial = obj.M;
                obj.M_angel = zeros(1, length(obj.M));
                if obj.change_type == 5 || obj.change_type == 6
                    [obj.x, obj.x_angel] = obj.ChangeMatrix(obj.x, obj.x_initial, obj.hub.lower, obj.hub.upper, obj.x_angel, 'x');
                    for i = 1:length(obj.M)
                        [obj.M{i}, obj.M_angel(i)] = obj.ChangeMatrix(obj.M{i}, obj.M_initial{i}, -1 * ones(1, obj.hub.D), 1 * ones(1, obj.hub.D), obj.M_angel(i), sprintf('M%d', i));
                        %obj.M{i} = orth(obj.M{i}); % todo É¾³ý
                        if ~all(round(obj.M{i} * obj.M{i}') == eye(size(obj.M{i})), 'all')
                            error('no orth matrix');
                        end
                    end
                end
                obj.x = obj.SetMinDistance(obj.x);
%                 obj.o = obj.x;
%                 obj.of = obj.hybrid_composition_func(obj.o, obj.gn_max, obj.sub_fun, obj.o, obj.sigma, obj.lambda, zeros(1, obj.gn_max), obj.M);
            end
           
            if obj.change_type == 7 || obj.change_type == 8
                obj.gn = randi(obj.proRand, obj.gn_max-1)+1;
                obj.selected_idx = [sort(randperm(obj.proRand, obj.gn_max, obj.gn)) obj.gn_max+1:obj.gn_max+obj.ln];
                if obj.change_type == 7 && obj.gn == 2
                    obj.is_add = true;
                end
            else
                obj.selected_idx = 1:obj.gn_max+obj.ln;
                obj.gn = obj.gn_max; 
            end
            
            obj.o = obj.x(obj.selected_idx(1:obj.gn), :);
            obj.of = obj.GetFits(obj.o);
%             obj.of = obj.of(obj.selected_idx(1:obj.gn));

            obj.his_o(end+1, :) = {obj.env, evaluated, obj.o};
            obj.his_of(end+1, :) = {obj.env, evaluated, obj.of};
        end
        
        function fits = GetFits(obj, decs)
            if obj.fun_num <= 4 % DF1
                dist = pdist2(decs, obj.x(obj.selected_idx, :));
                fits = max(obj.h(obj.selected_idx)' - obj.w(obj.selected_idx)' .* dist, [], 2);
            else % composition functions
%                 fits = obj.hybrid_composition_func(decs, obj.gn, obj.sub_fun, obj.o, obj.sigma, obj.lambda, zeros(1, obj.gn), obj.M);
                fits = obj.hybrid_composition_func(decs, obj.gn, obj.sub_fun(obj.selected_idx), obj.x(obj.selected_idx, :), obj.sigma(obj.selected_idx), obj.lambda(obj.selected_idx), zeros(1, obj.gn), obj.M(obj.selected_idx));
            end
        end
        
        function bestfit = GetBestFit(obj)
            bestfit = obj.his_of;
        end
        
        function bestpos = GetBestPos(obj)
            bestpos = obj.his_o;
        end
        
        function ChangeDynamic(obj)
%             disp(obj.env);
            if obj.hub.evaluated < obj.hub.evaluation
               %% change the environment
                obj.env = obj.env + 1;
                
                if obj.change_type == 7
                    if obj.is_add
                        obj.gn = obj.gn + 1;
                        if obj.gn == obj.gn_max
                            obj.is_add = false;
                        end
                    else
                        obj.gn = obj.gn - 1;
                        if obj.gn == 2
                            obj.is_add = true;
                        end
                    end
                    obj.selected_idx = 1:obj.gn;
                elseif obj.change_type == 8
                    obj.gn = randi(obj.proRand, obj.gn_max-1)+1;
                    obj.selected_idx = sort(randperm(obj.proRand, obj.gn_max, obj.gn));
                end
                
                if obj.change_type == 7 || obj.change_type == 8
                    obj.selected_idx = [obj.selected_idx obj.gn_max+1:obj.gn_max+obj.ln];
                    % adopt the C1 to change the parameters
                    type = obj.change_type;
                    obj.change_type = 1;
                    obj.env = obj.env - 1; % prevent from adding again
                    obj.ChangeDynamic();
                    obj.change_type = type;
                    return;
                end
                    
                    
                if obj.fun_num <= 4 % DF1
                    obj.h(obj.gn_max+1:end) = dynamic_change(obj.proRand, obj.h(obj.gn_max+1:end), obj.change_type, obj.h_bound(1), obj.h_bound(2), obj.h_s, obj.env);
                    obj.w = dynamic_change(obj.proRand, obj.w, obj.change_type, obj.w_bound(1), obj.w_bound(2), obj.w_s, obj.env);
                else % composition functions
                    for i = 1:length(obj.M)
                        [obj.M{i}, obj.M_angel(i)] = obj.ChangeMatrix(obj.M{i}, obj.M_initial{i}, -1 * ones(1, obj.hub.D), 1 * ones(1, obj.hub.D), obj.M_angel(i), sprintf('M%d', i));                        
                        %obj.M{i} = orth(obj.M{i}); % todo É¾³ý
                        if ~all(round(obj.M{i} * obj.M{i}') == eye(size(obj.M{i})), 'all')
                            error('no orth matrix');
                        end
                    end
                end
                
                % change the position of the optimal solutions
                [obj.x, obj.x_angel] = obj.ChangeMatrix(obj.x, obj.x_initial, obj.hub.lower, obj.hub.upper, obj.x_angel, 'x');
                obj.x = obj.SetMinDistance(obj.x);
                
%                 disp(obj.hub.evaluated);
                obj.o = obj.x(1:obj.gn, :);
                
%                 if obj.fun_num <= 4
%                     obj.of = obj.h(1:obj.gn);
%                 else
                    obj.of = obj.GetFits(obj.o);
%                 end
               %% save the changed postion and fitness
                obj.his_o(end+1, :) = {obj.env, obj.hub.evaluated, obj.o};
                obj.his_of(end+1, :) = {obj.env, obj.hub.evaluated, obj.of};
            end
        end
    end
    
    methods(Access = private)
        %% this function is used to calculate the fitness of combination function, such as CF1~CF4
        function res = hybrid_composition_func(obj, x, func_num, func, o, sigma, lambda, bias, M)
%             disp(obj.hub.evaluated);
            [ps,D] = size(x);
            weight = zeros(ps, length(func_num));
            for i = 1:func_num
                oo = repmat( o(i,:), ps, 1 );
                weight(:,i) = exp( -sum( (x-oo).^2, 2 )./2./( D * sigma(i)^2 ) );
            end

            max_weight = max(weight, [], 2);
            weight = (weight == max_weight) .* weight + (weight ~= max_weight) .* (weight .* (1-repmat(max_weight, 1, func_num).^10));

            if sum(weight,2) == 0
                weight = weight + 1;
            end
            
            weight = weight ./ repmat( sum(weight,2), 1, func_num );
            res = 0;
            for i = 1:func_num
                oo = repmat(o(i,:),ps,1);
                f = func{i}(((x-oo)./repmat(lambda(i,:),ps,1))*M{i});
%                 eval(['f = feval(func.f' int2str(i) ',((x-oo)./repmat(lambda(i,:),ps,1))*M.M' int2str(i) ');']);
                x1 = obj.hub.upper;
                f1 = func{i}((x1./lambda(i,:))*M{i});
%                 eval(['f1 = feval(func.f' int2str(i) ',(x1./lambda(i,:))*M.M' int2str(i) ');']);
                fit1 = 2000 .* f ./ f1;
                res = res + weight(:,i) .* ( fit1 + bias(i) );
            end
            res = -res;
        end
        
        %% this function is used to change the matrix after the environment change
        function [new_m, theta] = ChangeMatrix(obj, m, m_initial, lb, ub, old_theta, name)
% disp(obj.hub.evaluated);
            d = size(m, 2);
            if rem(d, 2) == 1
                l = d - 1;
            else
                l = d;
            end
            r = randperm(obj.proRand, d, l);
            if obj.change_type == 5 || obj.change_type == 6
                phi_min = 0;
                phi_max = pi/6;
                m = m_initial;
                if (obj.env >= 12)
                    r = obj.rs.(name)(rem(obj.env, 12)+1, :);
                end
            else
                phi_min = -pi;
                phi_max = pi;
            end
            if ~isfield(obj.rs, name)
                obj.rs.(name) = zeros(obj.maxEnv, l);
            end                
            obj.rs.(name)(obj.env+1, :) = r;
            if ismember(obj.change_type, [1,2,3])
                old_theta = 0;
            end
            theta = dynamic_change(obj.proRand, old_theta, obj.change_type, phi_min, phi_max, 1, obj.env);
            rotation_M = eye(d);
            for i = 1:l/2
                rotation_M(r(i*2-1),r(i*2-1)) = cos(theta);
                rotation_M(r(i*2),r(i*2)) = cos(theta);
                rotation_M(r(i*2-1),r(i*2)) = -sin(theta);
                rotation_M(r(i*2),r(i*2-1)) = sin(theta);
            end
            new_m  = m * rotation_M;
            new_m = boundary_check(new_m, lb, ub);
        end
        
        function pos = SetMinDistance(obj, pos)
            for i = 1:size(pos, 1)
                accept = false;
                
                if i == 1 
                    accept = true;
                else
                    dis = pdist2(pos(i, :), pos(1:i-1, :));
                    if all(dis(:) > obj.dpeaks) % accept the position
                        accept = true;
                    end
                end
                    
                while ~accept
                    v = rand(obj.proRand, 1, obj.hub.D) - 0.5;
                    len = sqrt(sum(v.^2, 2));
                    len = obj.dpeaks / len;
                    v = v * len;
                    pos(i,:) = pos(i, :) + v;
                    pos(i,:) = boundary_check(pos(i,:), obj.hub.lower, obj.hub.upper);

                    % check the distance
                    dis = pdist2(pos(i, :), pos(1:i-1, :));
                    if all(dis(:) > obj.dpeaks) % accept the position
                        accept = true;
                    end
                end
            end 
        end
    end
end