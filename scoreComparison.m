function results = scoreComparison(subNum)
% Get data (for one subject)
[modelFits, BICs] = comparePredictors(subNum);
numModels = size(modelFits, 3);

% This tells you which is the best model overall
bestModelInd = find(BICs == min(BICs));

% This compares each model for each vertex, within a given run,
% then selects the model chosen most often across all runs (per vertex).
% Best dimension-reduced to a set of ROIs somehow...
% so that you get a set of "ROI X is best explained by parameter Y"
[~, a] = max(modelFits, [], 3); % find max value of each model comparison
results = mode(a, 1); % which model wins most often per vertex?

% Now you have a 1xn vector of model numbers per vertex
% Subset them by ROI, and take the mode within the ROI?

% to do 

end