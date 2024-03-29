function ex3()
    debug = false;

    test_algorithms = {
        {@DynamicNichingDE, 3, 0}, ...
        {@DynamicNichingDE, 3, 1}, ...
        {@DynamicNichingDE, 3, 2}, ...
        {@DynamicNichingDE, 3, 3}, ...
    };

    test_problems = {
        {@DMMOP, 1, 1, 5}, ...
        {@DMMOP, 1, 2, 5}, ...
        {@DMMOP, 1, 3, 5}, ...
        {@DMMOP, 1, 4, 5}, ...
        {@DMMOP, 1, 5, 5}, ...
        {@DMMOP, 1, 6, 5}, ...
        {@DMMOP, 1, 7, 5}, ...
        {@DMMOP, 1, 8, 5}, ...
        {@DMMOP, 8, 1, 5}, ...
        {@DMMOP, 8, 2, 5}, ...
        {@DMMOP, 8, 3, 5}, ...
        {@DMMOP, 8, 4, 5}, ...
        {@DMMOP, 8, 5, 5}, ...
        {@DMMOP, 8, 6, 5}, ...
        {@DMMOP, 8, 7, 5}, ...
        {@DMMOP, 8, 8, 5}, ...
    };
    
    running(test_algorithms, test_problems, 3, debug);
%     % algorithm
%     for alg = 1:length(test_algorithms) 
%         for pro = 1:length(test_problems)
%             try
%                 delete(gcp('nocreate'));
%                 parpool('local',max_run);
%                 spmd(max_run)
%                     global_ = Global(test_algorithms{alg}, test_problems{pro}, 'run', labindex);
%                     result = global_.Start();
%                 end
%                 result = cat(1, result{1:end});
%                 pr = sum(result);
%                 pr = pr(1:end-1)/pr(end);
%                 dlmwrite(sprintf('./Data/EX3/A%d-P%d',  alg, pro), pr);
%             catch ErrorInfo
%                 pr = 1;
%                 dlmwrite(sprintf('./Data/EX3/A%d-P%d-ERRORERROR', alg, pro), pr);
%             end
%         end
%     end
end