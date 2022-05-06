function ex1()
    debug = false;

    test_algorithms = {
        {@DynamicNichingDE, 3, 1}, ...
        {@DynamicNichingPSO, 3, 1}, ...
        {@DynamicNichingCEP, 3, 1}, ...
        {@DynamicNichingGA, 3, 1}, ...
        {@DynamicNichingCSA, 3, 1}, ...
        {@DynamicNichingCMAES, 3, 1}, ...
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
        ...
        {@DMMOP, 1, 1, 5}, ...
        {@DMMOP, 1, 2, 5}, ...
        {@DMMOP, 1, 3, 5}, ...
        {@DMMOP, 1, 4, 5}, ...
        {@DMMOP, 1, 5, 5}, ...
        {@DMMOP, 1, 6, 5}, ...   
        {@DMMOP, 1, 7, 5}, ...   
        {@DMMOP, 1, 8, 5}, ...   
        ...
        {@DMMOP, 1, 1, 10}, ...
        {@DMMOP, 2, 1, 10}, ...
        {@DMMOP, 3, 1, 10}, ...
        {@DMMOP, 4, 1, 10}, ...
        {@DMMOP, 5, 1, 10}, ...
        {@DMMOP, 6, 1, 10}, ...
        {@DMMOP, 7, 1, 10}, ...
        {@DMMOP, 8, 1, 10}, ...     
    };
    
    running(test_algorithms, test_problems, 1, debug);
%     % algorithm
%     for alg = 1:length(test_algorithms) 
%         for pro = 5:length(test_problems)
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
%                 dlmwrite(sprintf('./Data/EX1/A%d-P%d',  alg, pro), pr);
%             catch ErrorInfo
%                 pr = 1;
%                 dlmwrite(sprintf('./Data/EX1/A%d-P%d-ERRORERROR', alg, pro), pr);
%             end
%         end
%     end
end