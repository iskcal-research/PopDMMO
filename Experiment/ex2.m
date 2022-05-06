function ex2()
    debug = false;

    test_algorithms = {
        {@DynamicNichingDE, 0, 1}, ...
        {@DynamicNichingDE, 1, 1}, ...
        {@DynamicNichingDE, 2, 1}, ...
        {@DynamicNichingDE, 3, 1}, ...
        {@DynamicNichingDE, 4, 1}, ...
    };

    test_problems = {
        {@DMMOP, 1, 1, 5}, ...
        {@DMMOP, 2, 1, 5}, ...
        {@DMMOP, 3, 1, 5}, ...
        {@DMMOP, 4, 1, 5}, ...
        {@DMMOP, 5, 1, 5}, ...
        {@DMMOP, 6, 1, 5}, ...
        {@DMMOP, 7, 1, 5}, ...
        {@DMMOP, 8, 1, 5}, ... 
    };
    
    running(test_algorithms, test_problems, 2, debug);
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
%                 dlmwrite(sprintf('./Data/EX2/A%d-P%d',  alg, pro), pr);
%             catch ErrorInfo
%                 pr = 1;
%                 dlmwrite(sprintf('./Data/EX2/A%d-P%d-ERRORERROR', alg, pro), pr);
%             end
%         end
%     end
end