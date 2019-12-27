classdef Individual < handle
    properties(SetAccess=private)
        dec;
        fit;
        con;
        add;
    end
    methods(Access = ?Hub)
       %% Constructor
        function obj = Individual(decs, fits, cons, adds)
            if nargin > 0 
                % each row represents a individual
                obj(size(decs, 1), 1) = Individual;
                % set these data to the individual arrays
                for i = 1:size(decs, 1)
                    obj(i).dec = decs(i, :);
                    obj(i).fit = fits(i, :);
                    obj(i).con = cons(i, :);
                end
                if ~isempty(adds)
                    for i = 1:size(decs, 1)
                        for prop = fieldnames(adds)
                            obj(i).add.(char(prop)) = adds.(char(prop))(i, :);
                        end
                    end
                end
            end
        end
    end
    
    methods        
        function value = decs(obj)
            value = cat(1, obj.dec);
        end
        
        function value = fits(obj)
            value = cat(1, obj.fit);
        end
        
        function value = cons(obj)
            value = cat(1, obj.con);
        end
        
        function varargout = adds(obj, props)
            if nargin == 1 % props is empty
                if ~isempty(obj(1).add)
                    props = fieldnames(obj(1).add);
                    data{:} = obj.adds(props);
                    pd = cell(1, 2*length(props));
                    pd(1:2:end) = props;
                    pd(2:2:end) = data;
                    varargout = struct(pd{:});
                else
                    varargout = cell(0);
                end
                return;
            end
            
            if ~isfield(obj(1).add, props)
                error("The addition information of individual does not contain the property: %s", props);
            end
            add_info = cat(1, obj.add);
            for i = 1:size(props, 1)
                varargout{i} = cat(1, add_info.(char(props(i))));
            end
        end
        
        function SetAdds(obj, adds)
            for i = 1:size(obj, 1)
                for prop = fieldnames(adds)
                    obj(i).add.(char(prop)) = adds.(char(prop))(i, :);
                end
            end
        end
    end
end