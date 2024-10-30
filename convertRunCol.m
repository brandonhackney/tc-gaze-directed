function [output] = convertRunCol(input)
% Input is an N x M predictor matrix with M predictors and N datapoints
% If the final column of input is a column of run numbers up to R,
% This function converts it to many binary columns.
% For example, if N = 2, it will turn [1;2] into [1 0;0 1]
% Also removes the second-to-last column, if it's a constant (i.e. all 1s)
% Output is then N * (M+R-2)
    
    if input(:,end-1) == ones(size(input,1),1)
        % Remove constant column to avoid near-singular matrices
        input(:,end-1) = [];
    end
    runNums = input(:,end);
    numRuns = max(runNums);
    runNums = repmat(runNums,1,numRuns);
    flag = [];
    for i = 1:numRuns
        runNums(runNums(:,i) ~= i,i) = 0;
        runNums(runNums(:,i) == i,i) = 1;
        if isequal(runNums(:,i), zeros(size(input,1),1))
            flag = [flag,i];
        end
    end
    % Remove columns of 0s (i.e. missing runs)
    runNums(:,flag) = [];
    
    output = [input(:,1:end-1) runNums];
end