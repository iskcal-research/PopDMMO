classdef (Abstract)Algorithm < handle
    properties(SetAccess=protected)
        algRand;        % the random generator
        hub;            % the hub object
        pop;            % the population used in the algorithm, often exposed by Record
    end
    methods(Access = public)
        function SetRand(obj, run)
            if isempty(obj.algRand)
                obj.algRand = RandStream.create('mt19937ar','seed',run);
            end
        end
    end
    methods(Access = public)             
        %% Constructor
        function obj = Algorithm()
            obj.hub = Hub.GetIns();
        end
        
        %% Used for the record in the global object
        function pop = Record(obj)
            pop = obj.pop;
        end
        
        %% Initialize the population
        function Initialize(obj)
            obj.pop = obj.hub.GetIndis(rand_indi(obj.algRand, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper));
        end
        
        %% Check whether the termination condition is satisfied
        function ter = Terminate(obj)
            if obj.hub.evaluated >= obj.hub.evaluation
                ter = true;
            else
                ter = false;
            end
        end
        
        %% Evolve the population in a generation
        function Evolve(obj)
            %% generate the offspring
            offspring = obj.GetOffspring(1:obj.hub.N);
            offspring = [offspring; obj.pop(length(offspring)+1:end)]; % the size may be smaller than pop
            %% select the next generation
            obj.pop = obj.Select(obj.pop, offspring);         
        end
        
        %% Generate the offspring individuals with the selected parents
        function offspring = GetOffspring(obj, select_idx)
            off_decs = obj.pop(select_idx).decs;
            offspring = obj.hub.GetIndis(off_decs);
        end
        
        %% Select the next population based on the parent population and its offspring
        function next_pop = Select(obj, parent, offspring)
            if length(offspring) < length(parent)
                offspring = [offspring; parent(length(offspring)+1:end)];
            end
            next_pop = parent;
            se_off = offspring.fits > parent.fits;
            next_pop(se_off) = offspring(se_off);
        end
        
        %% The whole process of the algorithms
        function Run(obj)
            obj.hub.Start();
        end
    end
end