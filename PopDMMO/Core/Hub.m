classdef Hub < handle
    properties(SetAccess = private)
        algorithm;              % the algoritm object
        problem;                % the problem object
        evaluated;              % the current consumed fittness evaluation
        run;                    % the current run index
        perform;                % the performance indicator
        perform_method;         % the performance calculator function
        debug;                  % the debug mode
        debug_pop;              % in debug mode, the population in each iteration would be saved
        command;                % perserve the parameter settings through input data
    end
    properties(SetAccess = ?Problem)
        encoding;               % the encoding of the problem
        D;                      % the dimension of the problem
        lower;                  % the lower boundary of the solution space
        upper;                  % the upper boundary of the solution space
        evaluation;             % the maximum number of evaluation
        freq;                   % the frequency of the dynamic problem, if the problem is not dynamic, freq is regarded as maximun evaluation.
    end
    properties(SetAccess = ?Algorithm)
        N;                      % the size of the population
        g;                      % the current generation (maybe write)
        result_pop;             % the key population used for calculate the performance (maybe write)
    end
    methods(Access = private)
        %% set the default values of the parameters by the hub objective
        function SetDefault(obj)
            obj.command = cell(0);
            obj.evaluated = 0;
            obj.evaluation = 15050;
            obj.freq = obj.evaluation;
            obj.run = 1;
            obj.g = 0;
            obj.result_pop = {};
            obj.perform = [];
            obj.perform_method = @get_peaks;
            obj.debug = false;
            obj.debug_pop = {};
            obj.encoding = "real";
            obj.D = 2;
            obj.lower = 0;
            obj.upper = 1;
            obj.N = 100;            
        end
        
        function SetCommandValue(obj, params)
            %% set the parameters by the command values
            propstr = {'evaluation', 'freq', 'run', 'perform_method', 'debug', 'D', 'lower', 'upper', 'N'};
            if ~isempty(params)
                check_props = cellfun(@(par) ischar(par) && ismember(par, propstr), params, 'UniformOutput', false);
                vector_props = cell2mat(check_props);
                for i = find(vector_props)
                    obj.command{end+1} = params{i};
                    obj.(params{i}) = params{i+1};
                end
            end
        end
        
        function SetRelated(obj)
            %% the random seed of the algorithm should be set
            obj.algorithm.SetRand(obj.run);
            
            %% after the parameters settings, the related paramters should be set
            % reset the boundary of the problem
            if obj.D > 1
%                 if isreal(obj.lower)
                if length(obj.lower) == 1
                    obj.lower = ones(1, obj.D) * obj.lower;
                elseif length(obj.lower) ~= obj.D
                    error("the lower boundary of the problem is illegal.");
                end
%                 if isreal(obj.upper)
                if length(obj.upper) == 1
                    obj.upper = ones(1, obj.D) * obj.upper;
                elseif length(obj.upper) ~= obj.D
                    error("the upper boundary of the problem is illegal.");
                end
            end
            % if the problem is not dynamic, freq is the same as evaluation
            if ~isa(obj.problem, 'DynamicProblem')
                obj.freq = obj.evaluation;
            end
        end
        
        function value = GetValue(obj, props, value)
            if ismember(props, obj.command)
                value = obj.(props);
            end
        end
    end
    
    methods
        %% Constructor
        function obj = Hub(varargin)
            % set the hub instance
            obj.GetIns(obj);
            % set the default values
            obj.SetDefault();
            % set the values by the parameters
            if nargin > 2
                %obj.SetCommandValue({varargin{3:end}});
                obj.SetCommandValue(varargin(3:end));
            end
            % set the problem object
            if length(varargin{2}) > 1
                obj.problem = varargin{2}{1}(varargin{2}{2:end});
            else
                obj.problem = varargin{2}();
            end
            % set the algorithm object
            if length(varargin{1}) > 1
                obj.algorithm = varargin{1}{1}(varargin{1}{2:end});
            else
                obj.algorithm = varargin{1}();
            end
            % set some parameters related to others
            obj.SetRelated();
            % initialize the problem
            obj.problem.Initialize(obj.evaluated);
        end
        
        %% Obtain the population by the two-dimension matrix
        function pop = GetIndis(obj, varargin)
            if isempty(obj.result_pop) || obj.g > obj.result_pop{end, 1}
                % check the rest evaluation number
                rest = obj.freq - rem(obj.evaluated, obj.freq);

                if nargin == 2
                    decs = varargin{1};
                    adds = [];
                elseif nargin == 3
                    decs = varargin{1};
                    adds = varargin{2};
                else
                    exit("call function GetIndividuals error");
                end
                n = min([rest, size(decs, 1)]);
                fits = obj.problem.GetFits(decs); % todo: calculate n individuals
                cons = obj.problem.GetCons(decs); % todo: constraint merge
                if ~isempty(adds)
                    for prop = fieldnames(adds)
                        if size(adds.(char(prop)), 1) > n
                            adds.(char(prop)) = adds.(char(prop))(1:n, :);
                        end
                    end
                end
                pop = Individual(decs(1:n, :), fits(1:n, :), cons(1:n, :), adds);
                obj.evaluated = obj.evaluated + n;
                if rem(obj.evaluated, obj.freq) == 0
                    obj.result_pop(end+1, :) = {obj.g, obj.evaluated, {}};
                    if isa(obj.problem, 'DynamicProblem') % the environment changes if the problem is dynamic
                        obj.problem.ChangeDynamic();
                    end
                end
            else
                pop = [];
            end
        end
        
        %% Start the standard process of the algorithm to solve the problem
        function perform = Start(obj)
            dynamic = isa(obj.algorithm, 'DynamicAlgorithm');
            obj.algorithm.Initialize();
            obj.g = obj.g + 1;
            obj.Debug(obj.algorithm.Record());
            while(~obj.algorithm.Terminate())
                obj.algorithm.Evolve();
                obj.Debug(obj.algorithm.Record());
                
                %% the final iteration in the static algorithm or in the current evironment of the dynamic algorithm
                if ~isempty(obj.result_pop) && isempty(obj.result_pop{end, 3})
                    assert(rem(obj.evaluated, obj.freq)==0, 'Current iteration is not the final iteration.');
                    obj.result_pop{end, 3} = obj.algorithm.Record();
                end
                
                obj.g = obj.g + 1;
%                 disp(obj.evaluated);
                
                if dynamic && obj.algorithm.CheckChange()
                    obj.algorithm.RespondChange();
                    obj.Debug(obj.algorithm.Record());
                end
            end
            perform = obj.perform_method(obj.result_pop, obj.problem.GetBestPos(), obj.problem.GetBestFit());
            obj.perform = perform;
        end
        
        %% Start the custom process of the algorithm to solve the problem
        function perform = StartByCustom(obj)
            obj.algorithm.Run();
            perform = obj.perform_method(obj.result_pop, obj.problem.GetBestFit(), obj.problem.GetBestPos());
            obj.perform = perform;
        end
        
        %% set the properties by other methods          
        function set.N(obj, value)
            obj.N = obj.GetValue('N', value);
        end
        
        function set.encoding(obj, value)
            obj.encoding = obj.GetValue('encoding', value);
        end
        
        function set.D(obj, value)
            obj.D = obj.GetValue('D', value);
        end
        
        function set.lower(obj, value)
            obj.lower = obj.GetValue('lower', value);
        end
        
        function set.upper(obj, value)
            obj.upper = obj.GetValue('upper', value);
        end
        
        function set.evaluation(obj, value)
            obj.evaluation = obj.GetValue('evaluation', value);
        end
        
        function set.freq(obj, value)
            obj.freq = obj.GetValue('freq', value);
        end
    end
    
    methods(Access = ?Algorithm)
       %% record data in iterations when record is set to true
        function Debug(obj, pop)
            if obj.debug
                obj.debug_pop(end+1, :) = {obj.g, obj.evaluated, pop};
            end
        end
    end
    
    methods(Static)
        %% the single hub object during the code executes
        function obj = GetIns(obj)
            persistent gIns;
            if nargin > 0
                gIns = obj;
            else
                obj = gIns;
            end
        end
    end
end