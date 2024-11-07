function varargout = splitRunRegressors(designMatrix, numPredictors)
% Assume a full design matrix with predictors of and not of interest,
% plus a trailing column indicating which run is which.
% Instead of having one column per source of nuisance for all runs,
% split into a separate column per run.
% So e.g. if you have two sources of nuisance, and three runs,
% you would end up with six new columns.

% Split the nuisance and run indicator columns off from the predictors
nuisance = designMatrix(:, numPredictors + 1 : end - 1);
runCol = designMatrix(:,end);
% Get some useful parameters
numRuns = numel(unique(runCol));
runMat = convertRunCol(runCol);
% Now use this to build out a really wide predictor matrix
newNuisance = [];
for r = 1:numRuns
    thisRunInd = runMat(:,r);
    thisRunDat = nuisance .* thisRunInd;
    newNuisance = [newNuisance, thisRunDat];
end % for each run
% Now deal out the exports
% But ignore the runMat here, because I already z-scored each run
if nargout > 0
    % Output 1 = predictors of interest
    varargout{1} = designMatrix(:, 1:numPredictors);
end
if nargout > 1
    % Output 2 = nuisance predictors
    varargout{2} = newNuisance;
end
if nargout > 2
    % Output 3 = run indicators
    varargout{3} = runMat;
end
end % function
    