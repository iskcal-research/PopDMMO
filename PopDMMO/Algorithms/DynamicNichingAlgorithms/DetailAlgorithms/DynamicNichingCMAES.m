classdef DynamicNichingCMAES < DynamicNichingBase
    properties(SetAccess=protected)
        %d = 0.1;
        %nc = 0.1;
        sub_len;
        sub_pop;
        % sub_size;
        eigen_g = 0;
        
        lambda;
        mu;
        xmean;
        zmean;
        sigma;
        weights;
        mueff;
        cc;
        cs;
        c1;
        cmu;
        damps;
        pc;
        ps;
        B;
        D;
        C;
        chiN;    
    end
    
    methods(Access=public)   
        function obj = DynamicNichingCMAES(nt, rt)
            obj = obj@DynamicNichingBase(nt, rt);
        end
        
        function Initialize(obj)
            obj.pop = obj.hub.GetIndis(rand_indi(obj.algRand, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper));
            
            obj.niching = obj.GetNiching();
            obj.InitializeParameter();
        end
        
        function Evolve(obj)
%             disp(obj.hub.evaluated);
%             obj.pop = [];
            start_pop = 1;
            for i = 1:obj.sub_len
%                 disp(obj.hub.evaluated);
                arz = randn(obj.algRand, obj.hub.D, obj.lambda{i});
                arx = zeros(obj.hub.D, obj.lambda{i});
                for k = 1:obj.lambda{i}
                    t = obj.xmean{i} + obj.sigma{i} * (obj.B{i} * obj.D{i} * arz(:, k));
                    arx(:, k) = (boundary_check(t', obj.hub.lower, obj.hub.upper))';
                    if ~all(arx(:, k)==t)
                        arz(:, k) = ((arx(:, k) - obj.xmean{i})./obj.sigma{i}) \ (obj.B{i} * obj.D{i});
                    end
                end
                
                obj.sub_pop{i} = obj.hub.GetIndis(arx');
%                 if ~isempty(obj.sub_pop{i})
%                     obj.pop(start_pop: start_pop+length(obj.sub_pop{i})-1) = obj.sub_pop{i};
%                 end
                
                if isempty(obj.sub_pop{i})
                    start_pop = start_pop + obj.lambda{i};
                    continue;
                end
                obj.pop(start_pop: start_pop+length(obj.sub_pop{i})-1) = obj.sub_pop{i};
                start_pop = start_pop + obj.lambda{i};
                

                [~, arindex] = sort(obj.sub_pop{i}.fits, 'descend');
                obj.xmean{i} = arx(:, arindex(1:obj.mu{i})) * obj.weights{i};
                obj.zmean{i} = arz(:, arindex(1:obj.mu{i})) * obj.weights{i};

                obj.ps{i} = (1-obj.cs{i})*obj.ps{i} + (sqrt(obj.cs{i}*(2-obj.cs{i})*obj.mueff{i})) * (obj.B{i} * obj.zmean{i});
                hsig = norm(obj.ps{i})/sqrt(1-(1-obj.cs{i})^(2*obj.hub.g))/obj.chiN{i} < 1.4+2/(obj.hub.D+1);
                obj.pc{i} = (1-obj.cc{i})*obj.pc{i} + hsig * sqrt(obj.cc{i}*(2-obj.cc{i})*obj.mueff{i}) * (obj.B{i}*obj.D{i}*obj.zmean{i});

                obj.C{i} = (1-obj.c1{i}-obj.cmu{i}) * obj.C{i} ... 
                    + obj.c1{i} * (obj.pc{i}*(obj.pc{i})' ... 
                    + (1-hsig) * obj.cc{i}*(2-obj.cc{i}) * obj.C{i}) ... 
                    + obj.cmu{i} ... 
                    * (obj.B{i}*obj.D{i}*arz(:,arindex(1:obj.mu{i}))) ...
                    * diag(obj.weights{i}) * (obj.B{i}*obj.D{i}*arz(:,arindex(1:obj.mu{i})))';

                obj.sigma{i} = obj.sigma{i} * exp((obj.cs{i}/obj.damps{i})*(norm(obj.ps{i})/obj.chiN{i} - 1)); 
                
                if obj.hub.g - obj.eigen_g > 1/(obj.c1{i}+obj.cmu{i})/obj.hub.D/10
                    obj.C{i} = triu(obj.C{i}) + triu(obj.C{i},1)';  
                    [obj.B{i},obj.D{i}] = eig(obj.C{i});           
                    obj.D{i} = diag(sqrt(diag(abs(obj.D{i}))));
                    obj.eigen_g = obj.hub.g;
                end
            end
            obj.hub.N = length(obj.pop);
        end
        
%         function niching = GetNiching(obj)
%             niching_tmp = GetNiching@DynamicNichingBase(obj);
%             locs = zeros(obj.hub.N, 1);
%             for i = 1:length(niching_tmp)
%                 for j = 1:length(niching_tmp{i})
%                     locs(niching_tmp{i}(j)) = i;
%                 end
%             end
%             selected = true(length(niching_tmp), 1);
%             for i = 1:length(niching_tmp)
%                 if length(niching_tmp{i}) == 1
%                     selected(i)=false;
%                     cur_dis = pdist2(obj.pop(niching_tmp{i}).dec, obj.pop.decs);
%                     cur_dis(niching_tmp{i}) = inf;
%                     [~, min_idx] = min(cur_dis);
%                     niching_tmp{locs(min_idx)}(end+1) = niching_tmp{i};
%                 end
%             end
%             niching_tmp = niching_tmp(selected);
%             niching = niching_tmp;
%         end
        
        function InitializeParameter(obj)
            %% initialize the parameters of each subpopulation
            obj.sub_len = length(obj.niching);
            obj.sub_pop = cell(0);
            
            % obj.sub_size = zeros(1, obj.sub_len);
            obj.lambda = cell(1, obj.sub_len);
            obj.mu = cell(1, obj.sub_len);
            obj.sigma = cell(1, obj.sub_len);
            obj.weights = cell(1, obj.sub_len);
            obj.mueff = cell(1, obj.sub_len);
            obj.cc = cell(1, obj.sub_len);
            obj.cs = cell(1, obj.sub_len);
            obj.c1 = cell(1, obj.sub_len);
            obj.cmu = cell(1, obj.sub_len);
            obj.damps = cell(1, obj.sub_len);
            obj.pc = cell(1, obj.sub_len);
            obj.ps = cell(1, obj.sub_len);
            obj.B = cell(1, obj.sub_len);
            obj.D = cell(1, obj.sub_len);
            obj.C = cell(1, obj.sub_len);
            obj.chiN = cell(1, obj.sub_len);  
            obj.xmean = cell(1, obj.sub_len);
            obj.zmean = cell(1, obj.sub_len);
            
            dim = obj.hub.D;
            for i = 1:obj.sub_len
                obj.lambda{i} = 4 + floor(3*log(dim));
                obj.mu{i} = obj.lambda{i} / 2;
                obj.sigma{i} = 0.5;
                obj.weights{i} = log(obj.mu{i}+1/2)-log(1:obj.mu{i})';
                obj.mu{i} = floor(obj.mu{i});
                obj.weights{i} = obj.weights{i} / sum(obj.weights{i});
                obj.mueff{i} = sum(obj.weights{i})^2/sum(obj.weights{i}.^2);

                obj.cc{i} = (4+ obj.mueff{i}/dim) / (dim+4 + 2* obj.mueff{i}/dim);
                obj.cs{i} = (obj.mueff{i} + 2)/(dim + obj.mueff{i} + 5);
                obj.c1{i} = 2 / ((dim+1.3)^2 + obj.mueff{i}); 
                obj.cmu{i} = 2 * (obj.mueff{i} - 2 + 1 / obj.mueff{i}) / ((dim + 2) ^2 + 2 * obj.mueff{i}/2);
                obj.damps{i} = 1 + 2 * max(0, sqrt((obj.mueff{i} - 1)/(dim + 1)) - 1) + obj.cs{i};

                obj.pc{i} = zeros(dim,1);
                obj.ps{i} = zeros(dim,1); 
                obj.B{i} = eye(dim); 
                obj.D{i} = eye(dim);
                obj.C{i} = obj.B{i} * obj.D{i} * (obj.B{i} * obj.D{i})';
                obj.chiN{i} = dim^0.5 * (1 - 1 / (4 * dim) + 1 / (21 * dim^2));
                
                obj.xmean{i} = (obj.pop(obj.niching{i}(1)).dec)';
%                 sub_pop = obj.pop(obj.niching{i}(1:obj.mu{i})); 
%                 [~, sort_idx] = sort(sub_pop.fits, 'descend');
%                 obj.xmean{i} = mean(sub_pop(sort_idx(1:obj.mu{i})).decs, 1)';
%                 obj.zmean{i} = zeros(obj.hub.D, 1);
            end 
        end
        
        function RespondChange(obj)
            RespondChange@DynamicNichingBase(obj);
            obj.niching = obj.GetNiching();
            obj.InitializeParameter();
            obj.eigen_g = obj.hub.g+1;
        end
    end
end