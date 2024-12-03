function plotPredCorr(fullPred, subNum, predLabels)
% Given the "full" predictor matrix (all runs) for one subject,
% which comes from comparePredictors(),
% plot the correlation matrix for the predictors of interest.
% (Though apparently everyone had the same stimulus order,
% so they should all be exactly the same.)

% Define the expected list of predictors. (probably should be an input).
% predLabels = {'Motion', 'Interactivity', 'Top-Down Gaze', 'Understanding', 'StimTiming'};
numPred = length(predLabels);
assert(width(fullPred) >= numPred, 'Not enough columns in this predictor matrix!');

figure();
h = heatmap(corr(fullPred(:,1:numPred), 'rows', 'complete'));
h.XDisplayLabels = predLabels;
h.YDisplayLabels = predLabels;

subID = sprintf('sub-%02.f', subNum);
title(subID);

end