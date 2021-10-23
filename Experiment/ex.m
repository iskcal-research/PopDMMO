function ex()
    func = { ...
%         @function_plot, ...
         @ex1, ...
        @ex2, ...
        @ex3, ...
    };

    for i = 1:length(func)
        clearvars -except func i;
        clc;
        fprintf('Experiment (%s) is running\n', func2str(func{i}));
        func{i}();
    end
end