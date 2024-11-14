function labelStruct = scoreComparison(modelFits, labelStruct, BICs)
% Run from within comparePredictors
% Take the smaller stack of model fits (i.e. one per ROI, not whole-brain)
% Return the best model (and associated R2) per ROI
numModels = size(modelFits, 3);

% This tells you which is the best model overall
bestModelInd = find(BICs == min(BICs));

% This compares each model for each ROI, within a given run,
% then selects the model chosen most often across all runs (per ROI).
% so that you get a set of "ROI X is best explained by parameter Y"
[~, a] = max(modelFits, [], 3); % find max value of each model comparison
results = mode(a, 1); % which model wins most often per ROI?
clear a

% Now you have a 1xn vector of model numbers per ROI
% Write that, and the associated R2 (or average R2 I guess), to the struct
assert(length(labelStruct) == length(results), 'Model fits and labels have unequal number of ROIs!');
roimeans = [];
roimaxes = [];
for i = 1:length(results) 
    labelStruct(i).BestModel = results(i);
    x = modelFits(:, i, results(i));
    labelStruct(i).R2 = mean(x, 1); % also export R2 val
    % Above, we're taking the best model per run per ROI,
    % then averaging the max R2 across all runs of that model.
    % (As opposed to subsetting to just the specific runs when that model
    % was the best).
    
    % Temp: report
    fprintf(1, '%s:\tMax R\x00B2 value = %0.4f\tMean R\x00B2 value = %0.4f\n', labelStruct(i).Label, max(x), mean(x));
    roimeans(i) = mean(x);
    roimaxes(i) = max(x);
end

fprintf(1, '\n\nMean Max = %0.4f\tMax Mean = %0.4f\n', mean(roimaxes, 'omitnan'), max(roimeans));
fprintf(1, 'Max Max = %0.4f\tMean Mean = %0.4f\n', max(roimaxes), mean(roimeans, 'omitnan'));