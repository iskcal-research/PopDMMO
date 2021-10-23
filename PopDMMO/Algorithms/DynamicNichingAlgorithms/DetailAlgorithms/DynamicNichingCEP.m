classdef DynamicNichingCEP < DynamicNichingBase
   properties(SetAccess=protected)
        tau;
        tau1;
        q = 10;
    end
    methods
        function obj = DynamicNichingCEP(nt, rt)
            obj = obj@DynamicNichingBase(nt, rt);
        end
        
        function Initialize(obj)
            obj.tau = 1/(sqrt(2*sqrt(obj.hub.D)));
            obj.tau1 = 1/(sqrt(2*obj.hub.D));
            obj.pop = obj.hub.GetIndis(rand_indi(obj.algRand, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper), struct("eta", 3 * ones(obj.hub.N, obj.hub.D)));
        end
        
        function offspring = GetOffspring(obj, select_idx)
            pop_decs = obj.pop(select_idx).decs;
            eta = obj.pop(select_idx).adds("eta");
            
            off_decs = pop_decs + eta .* randn(obj.algRand, size(pop_decs));
            off_decs = boundary_check(off_decs, obj.hub.lower, obj.hub.upper);
            eta_new = eta .* exp(obj.tau1 .* repmat(randn(obj.algRand, size(pop_decs, 1), 1), 1, obj.hub.D) +  obj.tau .* randn(obj.algRand, size(pop_decs)));
            offspring = obj.hub.GetIndis(off_decs, struct("eta", eta_new));
        end
        
        function next_pop = Select(obj, pop, offspring)
            indis = [pop; offspring];
            idx = tourment_select(obj.algRand, indis.fits, length(pop), obj.q);
            next_pop = indis(idx);
        end
        
         %% dynamic response 
        function RespondChange(obj)
            RespondChange@DynamicNichingBase(obj);
            obj.pop.SetAdds(struct("eta", 3 * ones(obj.hub.N, obj.hub.D)));
        end
    end
end