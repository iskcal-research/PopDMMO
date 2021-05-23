function running(test_algorithms, test_problems, num)
    max_run = 30;
    
%     test_algorithms = {{@DynamicNichingDE, 1, 1}};
%     test_problems = {{@DMMOP, 1, 1, 5}};
%     num = 2;

    if ~exist(sprintf('./Data/EX%d', num),'dir')
        mkdir(sprintf('./Data/EX%d', num));
    end

    % algorithm
    for alg = 1:length(test_algorithms)
        for pro = 1:length(test_problems)
            fprintf('Algorithm (%d) to solve Problem (%d) is running\n', alg, pro);
            try
                delete(gcp('nocreate'));
                parpool('local',max_run);
                spmd(max_run)
                    hub = Hub(test_algorithms{alg}, test_problems{pro}, 'run', labindex);
                    result = hub.Start();
                end
                result = cat(1, result{1:end});
                pr = sum(result);
                pr = pr(1:end-1)/pr(end);
                dlmwrite(sprintf('./Data/EX%d/A%d-P%d', num, alg, pro), pr);
            catch exception
                fileID = fopen(sprintf('./Data/EX%d/A%d-P%d-ERRORERROR', num, alg, pro), 'w');
                fprintf(fileID, sprintf('The Algorithm (%d) to solve Problem (%d) in Experiment (%d) is error.\r\n', alg, pro, num));
                fprintf(fileID, [exception.identifier '\r\n']);
                fprintf(fileID, [exception.message '\r\n']);
                fprintf(fileID, 'ErrorCause:\r\n');
                                
                fprintf(fileID, 'StackInfo:\r\n');
                for i = 1:length(exception.stack)
                    fprintf(fileID, sprintf('\t Line (%d), Function (%s), File (%s)\r\n', exception.stack(i).line, exception.stack(i).name, strrep(exception.stack(i).file, '\', '\\')));
                end
                
                for i = 1: length(exception.cause)
                    inner_exception = exception.cause{i};
                    fprintf(fileID, ['\t' inner_exception.identifier , '----', inner_exception.message '\r\n']);
                    for j = 1:length(inner_exception.stack)
                        fprintf(fileID, sprintf('\t\t Line (%d), Function (%s), File (%s)\r\n', inner_exception.stack(j).line, inner_exception.stack(j).name, strrep(inner_exception.stack(j).file, '\', '\\')));
                    end 
                end

                
                fclose(fileID);
%                 pr = 1;
%                 dlmwrite(sprintf('./Data/EX%d/A%d-P%d-ERRORERROR', num, alg, pro), pr);
            end
        end
    end
end