function [fullPred, predList] = getPredictorStack(subNum, predList)

numRuns = countRuns(subNum, 'tricopa');
fullPred = [];

% Now go:
fprintf(1, 'Subject %i: ', subNum);
fprintf(1, 'Getting predictors for %i runs\n', numRuns);

for r = 1:numRuns
    % Get the design matrix for each run
    [tmpPred, predList] = getSDM(subNum, r, predList);
    tmpPred = zscore(tmpPred);
    % Insert an intercept column
    tmpPred = [ones(height(tmpPred), 1), tmpPred];
    % Insert run num as a trailing column, for later subsetting
    tmpPred(:, end+1) = r;
    fullPred = [fullPred; tmpPred]; % stack runs
end
fprintf(1, 'Done getting predictors for subject %i.\n', subNum);
end