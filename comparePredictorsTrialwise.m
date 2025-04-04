function [results, roiLabels, predList] = comparePredictorsTrialwise(varargin)
% 1. Aggregate MRI data for all runs of a specific task 
% 2. Regress out effects of nuisance parameters (like CSF signal)
% 3. Generate a predictor matrix treating each trial as a unique condition
% 4. Estimate effect of trial timing, one against all, on residuals
% 5. Extract the beta for each trial, and the average parameter values
% 6. Correlate the beta series with the parameter series
% 7. See which parameter has the highest r score per ROI

numSubs = 0;
if nargin > 0
    % Allow a list of subjects to iterate over (even just one)
    subList = varargin{1};
    numSubs = length(subList);
end
if numSubs == 0
    % if no input OR if empty input
    numSubs = 30;
    subList = 1:numSubs;
end

hem = 2; % for RH
for subNum = subList
    [dataStack,  roiLabels, fullPred, predList] = getDataStack(subNum, hem);
    dataStack = stackRuns(dataStack);
    trialPred = getTrialTiming(subNum);
    
    fullPred(:,[1, end]) = []; % drop the intercept and run number cols
    
    % Now for the meat of it.
    % trialPred is n * 2 * 100, so you need to:
    % - loop over all 100 trials
    % - estimate the GLM for that page only
    % - extract the beta for the on-trial only (i.e. index 1)
    % - store that in a vector for all trials
    % - get a parameter-average value for each predictor
    % - store those in a matrix for all trials
    % then repeat for all subjects, so make sure all outputs are indexed.
    % In the end, find the correlation between the betas and parameters.
    numTrials = size(trialPred, 3); % should be 100
    betaSeries = [];
    paramSeries = [];
    for t = 1:numTrials
        betas = simpleGLM(dataStack, trialPred(:,:,t));
        betaSeries(t,:) = betas(2,:);
        % Reduce parameters down to a single number
        y = fullPred .* trialPred(:,1,t);
        y = y(trialPred(:,1,t) ~= 0,:);
        paramSeries(t,:) = mean(y, 1);
    end
    % Hopefully this works out
    results(subNum,:,:) = corr(betaSeries, paramSeries);
end

% So now we've got a 3D matrix of i subjects * j ROIs * k predictors
% That is the format expected by scoreComparison,
% except it expects R^2s (not r), and an extra page of "full model" data.
% Sending it what we have here won't make it crash,
% but some of the output will be wrong and/or meaningless.
% May want to write a separate function specifically for trial-based data?
% Or write in an exception that prints things differently etc.
scoreComparison(results, roiLabels, predList);
boxplotROIs(results, roiLabels, predList);
end