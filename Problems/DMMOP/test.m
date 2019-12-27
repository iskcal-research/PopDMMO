classdef test < DynamicProblem
    properties(SetAccess=private)
        fun_num = 1;        % function number, if it less than 5, the problem is based on MPB, otherwise, it is based on MMPOP.
        change_type = 1;
        o;                  % the best position
        of;                  % the best fitness
        gn = 4;
        
%        dpeak = 0.1;
        % the following four parameters are used in MPB
        ln = 0;
        x;
        h;
        w;
        h_s = 7.0;
        w_s = 1.0;
        h_bound = [30, 70];
        w_bound = [1, 12];
        x_initial;
        
        % the following one parameter is used in MMOP
        M;
        M_initial;
        o_initial;
        sub_fun;
        lambda;
        sigma;
        
        % history best position and fitness;
        his_o = {};
        his_of = {};
    end
    methods       
        function obj = DMMOP(fun_num, change_type, D, varargin)
            obj.fun_num = fun_num;
            obj.change_type = change_type;
            obj.global_.D = D;
            obj.global_.lower = -5 .* ones(1, D);
            obj.global_.upper = 5 .* ones(1, D);
            obj.global_.freq = 3000 * D;
            obj.global_.evaluation = 60 * obj.global_.freq;
            x_range = 10;
            
            %% set the parameters in DMMOP
            if fun_num <= 4 % DMMOP is based on MPB
                if nargin > 3 % set the number of the global peaks
                    obj.gn = varargin{1};
                end
                if nargin > 4 % set the number of the local peaks
                    obj.ln = varargin{2};
                end
                
                %% set parameters of MPB
                h_global = 75;
                switch(fun_num)
                    case 1 % absolutely random 
                        obj.h = [ones(obj.gn, 1) .* h_global; rand(obj.proRand, obj.ln, 1) .* (obj.h_bound(2) - obj.h_bound(1)) + obj.h_bound(1)];
                        obj.w = rand(obj.proRand, obj.gn+obj.ln, 1) .* (obj.w_bound(2) - obj.w_bound(1)) + obj.w_bound(1);
                        obj.x = rand(obj.proRand, obj.gn+obj.ln, obj.global_.D) .* (obj.global_.upper - obj.global_.lower) + obj.global_.lower;
                    case 2 % simulate Shubert function
                        x_dis = x_range * 0.3;
                        if rem(obj.gn, 2) == 1
                            obj.gn = obj.gn + 1;
                        end
                        if obj.global_.D > 1
                            obj.ln = obj.gn;
                        else
                            obj.ln = 0;
                        end
                        obj.h = [ones(obj.gn, 1) .* h_global; rand(obj.proRand, obj.ln, 1) .* (obj.h_bound(2) - obj.h_bound(1)) + obj.h_bound(1)];
                        obj.w = ones(obj.gn+obj.ln, 1) .* 5;
                        x_original = rand(obj.proRand, 1, obj.global_.D) .* ((obj.global_.upper-x_dis) - (obj.global_.lower+x_dis)) + (obj.global_.lower+x_dis);
                        if obj.gn > 2*(1+2^obj.global_.D)
                            error("too many global peaks in function 2");
                        end
                        for i = 1:obj.global_.D
                            x_original(end+1, :) = x_original(1, :);
                            x_original(end, i) = x_original(end, i) - x_dis;
                            x_original(end+1, :) = x_original(1, :);
                            x_original(end, i) = x_original(end, i) + x_dis;
                        end
                        x_original = x_original(1:obj.gn/2, :);
                        dir = repmat([-1;1] .* 0.05 .* x_range, obj.gn/2, obj.global_.D);
                        obj.x = repelem(x_original, 2, 1) + dir;
                        if obj.global_.D > 1
                            dir(:, 2) = -dir(:, 2);
                            obj.x = [obj.x; repelem(x_original, 2, 1) + dir];
                        end
                    case 3 % simulate Vincent function
                        obj.ln = 0;
                        
                        obj.h = ones(obj.gn, 1) .* h_global;
                        obj.w = ones(obj.gn+obj.ln, 1) .* 5;
                        x_original = rand(obj.proRand, 1, obj.global_.D) .* ((obj.global_.upper-x_range * 0.2) - (obj.global_.lower+x_range * 0.2)) + (obj.global_.lower+x_range * 0.2);
                        if obj.gn > 1+2^obj.global_.D
                            error("too many global peaks in function 3");
                        elseif obj.gn < 3
                            error("too less global peaks in function 3");
                        end
                        for i = 1:obj.global_.D
                            x_original(end+1, :) = x_original(1, :);
                            x_original(end, i) = x_original(end, i) - x_range * 0.15;
                            obj.w = [obj.w; 2];
                            x_original(end+1, :) = x_original(1, :);
                            x_original(end, i) = x_original(end, i) + x_range * 0.3;
                            obj.w = [obj.w; 9];
                        end
                        obj.x = x_original(1:obj.gn, :);
                        obj.w = obj.w(1:obj.gn, :);
                    case 4 % simulate modified Rastrigin function
                        obj.ln = 0;
                        obj.h = ones(obj.gn, 1) .* h_global;
                        obj.w = ones(obj.gn+obj.ln, 1) .* 5;
                        x_dis = 1;
                        x_original = rand(obj.proRand, 1, obj.global_.D) .* ((obj.global_.upper-x_dis) - (obj.global_.lower+x_dis)) + (obj.global_.lower+x_dis);
                        if obj.gn > 1+2^obj.global_.D
                            error("too many global peaks in function 4");
                        end
                        for i = 1:obj.global_.D
                            x_original(end+1, :) = x_original(1, :);
                            x_original(end, i) = x_original(end, i) - x_dis;
                            x_original(end+1, :) = x_original(1, :);
                            x_original(end, i) = x_original(end, i) + x_dis;
                        end
                        obj.x = x_original(1:obj.gn, :);
                    otherwise
                        error('fun_num in DMMOP is illegal.');
                end
                
                obj.x_initial = obj.x;
                if obj.change_type == 5 || obj.change_type == 6
                    obj.h = [ones(obj.gn, 1); dynamic_change(obj.proRand, obj.h(obj.gn+1:end), obj.change_type, obj.h_bound(1), obj.h_bound(2), obj.h_s, obj.env)];
                    obj.w = dynamic_change(obj.proRand, obj.w, obj.change_type, obj.w_bound(1), obj.w_bound(2), obj.w_s, obj.env);                    
                    obj.x = obj.ChangeMatrix(obj.x, obj.x, obj.global_.lower, obj.global_.upper);
                end
                obj.of = obj.h(1:obj.gn);
                obj.o = obj.x(1:obj.gn, :);
            else % DMMOP is based on MMOP
                switch(fun_num)
                    case 5  % CF1
                        obj.sub_fun = {@FGrienwank, @FGrienwank, @FWeierstrass, @FWeierstrass, @FSphere, @FSphere};
                        obj.gn = length(obj.sub_fun);
                        obj.sigma = ones(1, obj.gn);
                        obj.lambda = [1, 1, 8, 8, 1/5, 1/5]';
                        obj.M = cell(obj.gn, 1);
                        obj.M(1:end) = {diag(ones(1, obj.global_.D))};
                    case 6 % CF2
                        obj.sub_fun = {@FRastrigin, @FRastrigin, @FWeierstrass, @FWeierstrass, @FGrienwank, @FGrienwank, @FSphere, @FSphere};
                        obj.gn = length(obj.sub_fun);
                        obj.sigma = ones(1, obj.gn);
                        obj.lambda = [1, 1, 10, 10, 1/10, 1/10, 1/7, 1/7]';
                        obj.M = cell(obj.gn, 1);
                        obj.M(1:end) = {diag(ones(1, obj.global_.D))};
                    case 7 %CF3
                        obj.sub_fun = {@FEF8F2, @FEF8F2, @FWeierstrass, @FWeierstrass, @FGrienwank, @FGrienwank};
                        obj.gn = length(obj.sub_fun);
                        obj.sigma = [1,1,2,2,2,2];
                        obj.lambda = [1/4, 1/10, 2, 1, 2, 5]';
                        if ~ismember(obj.global_.D, [2, 3, 5, 10, 20])
                            error('the dimension is illegal in function 7.');
                        end
                        file = fullfile(fileparts(mfilename('fullpath')), sprintf('./Data/CF3_M_D%d', obj.global_.D));
                        t = load(file);
                        obj.M = t.M;
                    case 8 %CF4
                        obj.sub_fun = {@FRastrigin, @FRastrigin, @FEF8F2, @FEF8F2, @FWeierstrass, @FWeierstrass, @FGrienwank, @FGrienwank};
                        obj.gn = length(obj.sub_fun);
                        obj.sigma = [1,1,1,1,1,2,2,2];
                        obj.lambda = [4, 1, 4, 1, 1/10, 1/5, 1/10, 1/40]';
                        if ~ismember(obj.global_.D, [2, 3, 5, 10, 20])
                            error('the dimension is illegal in function 8.');
                        end
                        file = fullfile(fileparts(mfilename('fullpath')), sprintf('./Data/CF4_M_D%d', obj.global_.D));
                        t = load(file);
                        obj.M = t.M;
                    otherwise
                        error('fun_num in DMMOP is illegal.');
                end
                file = fullfile(fileparts(mfilename('fullpath')),'./Data/optima');
                t = load(file);
                
                obj.o = t.o(1:obj.gn, 1:obj.global_.D);
                obj.o_initial = obj.o;
                obj.M_initial = obj.M;
                if obj.change_type == 5 || obj.change_type == 6
                    obj.o = obj.ChangeMatrix(obj.o, obj.o_initial, obj.global_.lower, obj.global_.upper);
                    for i = 1:length(obj.M)
                        obj.M{i} = obj.ChangeMatrix(obj.M{i}, obj.M_initial{i}, -1 * ones(1, obj.global_.D), 1 * ones(1, obj.global_.D));
                        obj.M{i} = orth(obj.M{i});
                    end
                end
                
                obj.of = obj.hybrid_composition_func(obj.o, obj.gn, obj.sub_fun, obj.o, obj.sigma, obj.lambda, zeros(1, obj.gn), obj.M);
            end
            
            obj.his_o(end+1, :) = {obj.env, obj.global_.evaluated, obj.o};
            obj.his_of(end+1, :) = {obj.env, obj.global_.evaluated, obj.of(1)};
        end
        
        function fits = GetFits(obj, decs)
            if obj.fun_num <= 4 % MBP
                dist = pdist2(decs, obj.x);
                fits = max(obj.h' - obj.w' .* dist, [], 2);
            else % MMOP
                fits = obj.hybrid_composition_func(decs, obj.gn, obj.sub_fun, obj.o, obj.sigma, obj.lambda, zeros(1, obj.gn), obj.M);
            end
        end
        
        function bestfit = GetBestFit(obj)
            bestfit = obj.his_of;
        end
        
        function bestpos = GetBestPos(obj)
            bestpos = obj.his_o;
        end
        
        function DynamicChange(obj)
            if obj.global_.evaluated < obj.global_.evaluation
               %% ochange the environment
                obj.env = obj.env + 1;
                if obj.fun_num <= 4 % MPB
                    obj.h(obj.gn+1:end) = dynamic_change(obj.proRand, obj.h(obj.gn+1:end), obj.change_type, obj.h_bound(1), obj.h_bound(2), obj.h_s, obj.env);
                    obj.w = dynamic_change(obj.proRand, obj.w, obj.change_type, obj.w_bound(1), obj.w_bound(2), obj.w_s, obj.env);
                  %% change the position of the optimal solutions
                    obj.x = obj.ChangeMatrix(obj.x, obj.x_initial, obj.global_.lower, obj.global_.upper);
                    
                    obj.o = obj.x(1:obj.gn, :);
                    obj.of = obj.h(1:obj.gn);
                else % MMOP
                    obj.x = obj.ChangeMatrix(obj.o, obj.o_initial, obj.global_.lower, obj.global_.upper);
                    for i = 1:length(obj.M)
                        obj.M{i} = obj.ChangeMatrix(obj.M{i}, obj.M_initial{i}, -1 * ones(1, obj.global_.D), 1 * ones(1, obj.global_.D));
                        obj.M{i} = orth(obj.M{i});
                    end
                end

               %% save the changed postion and fitness
                obj.his_o(end+1, :) = {obj.env, obj.global_.evaluated, obj.o};
                obj.his_of(end+1, :) = {obj.env, obj.global_.evaluated, obj.of(1)};
            end
        end
    end
    
    methods(Access = private)
        %% this function is used to calculate the fitness of combination function, such as CF1~CF4
        function res = hybrid_composition_func(obj, x, func_num, func, o, sigma, lambda, bias, M)
%             disp(obj.global_.evaluated);
            [ps,D] = size(x);
            weight = zeros(ps, D);
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
                x1 = obj.global_.upper;
                f1 = func{i}((x1./lambda(i,:))*M{i});
%                 eval(['f1 = feval(func.f' int2str(i) ',(x1./lambda(i,:))*M.M' int2str(i) ');']);
                fit1 = 2000 .* f ./ f1;
                res = res + weight(:,i) .* ( fit1 + bias(i) );
            end
            res = -res;
        end
        
        %% this function is used to change the matrix after the environment change
        function new_m = ChangeMatrix(obj, m, m_initial, lb, ub)
            if obj.change_type == 4
                new_m = dynamic_change(obj.proRand, m, obj.change_type, lb, ub, [], obj.env);
            else
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
                else
                    phi_min = -pi;
                    phi_max = pi;
                end
                theta = dynamic_change(obj.proRand, 0, obj.change_type, phi_min, phi_max, 1, obj.env);
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
        end
    end
end