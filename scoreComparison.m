function labelStruct = scoreComparison(modelFits, labelStruct, predictorList)
% Run from within comparePredictors
% Take the smaller stack of model fits (i.e. one per ROI, not whole-brain)
% Return the best model (and associated R2) per ROI
numModels = size(modelFits, 3);
numPreds = length(predictorList);

% This compares each model for each ROI, within a given run,
% then selects the model chosen most often across all runs (per ROI).
% so that you get a set of "ROI X is best explained by parameter Y"
[~, a] = max(modelFits, [], 3); % find max value of each model comparison
results = mode(a, 1); % across runs, which model wins most often per ROI?
clear a

% An alternative method is to average the results across subjects,
% then decide which model has the maximum value per ROI.
modelAvgs = mean(modelFits, 1);

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
    [~, bestModel(i)] = max(modelAvgs(:,i,:));
    [~, worstModel(i)] = min(modelAvgs(:,i,:));
    % Temp: report
    fprintf(1, '%s:\t',labelStruct(i).Label);
%     fprintf(1, '\tMax R\x00B2 value = %0.3f,', max(x));
%     fprintf(1, '\tMean R\x00B2 value = %0.3f,', mean(x));
%     fprintf(1, '\tBest Model = %i', bestModel(i));
%     fprintf(1, '\tWorst Model = %i\n', worstModel(i));
    
    % Variance partitioning:
    % R2 of full model MINUS R2 of model missing one predictor
    % Significance level comes from bootstrapping 10,000 shuffled iterations
    fprintf(1, 'Variance attributable to:');
    for p = 1:numPreds+1
        if p == numPreds + 1
            fprintf('\t%s: %0.4f', 'FullModel', modelAvgs(:,i,end));
%             fprintf('\t%s: %0.4f', 'FullModel', 1); % meaningless if x-y/x
        else
            predName = predictorList{p};
            fullmodel = modelAvgs(:,i,end);
            limmodel = modelAvgs(:,i,p);
%             fprintf('\t%s: %0.4f', predName, fullmodel - limmodel);
            fprintf('\t%s: %0.4f', predName, (fullmodel - limmodel) / fullmodel);
        end
    end
    fprintf(1, '\n');
    roimeans(i) = mean(x);
    roimaxes(i) = max(x);
end

fprintf(1, '\n\nMean Max = %0.4f\tMax Mean = %0.4f\n', mean(roimaxes, 'omitnan'), max(roimeans));
fprintf(1, 'Max Max = %0.4f\tMean Mean = %0.4f\n', max(roimaxes), mean(roimeans, 'omitnan'));

% Most consistent ROIs (i.e. peak average across all subjects)
[v, i] = maxk(roimeans, 3);
n = {labelStruct(i).Label};
fprintf(1, '\nThe three most consistent ROIs are:\n');
fprintf(1, '\t%s (mean R2 = %0.4f),\n', n{1}, v(1));
fprintf(1, '\t%s (mean R2 = %0.4f), and \n', n{2}, v(2));
fprintf(1, '\t%s (mean R2 = %0.4f).\n', n{3}, v(3));

% Single-subject best model fit
[v, i] = maxk(roimaxes, 3);
n = {labelStruct(i).Label};
fprintf(1, '\nThe three highest single-subject fits are in:\n');
fprintf(1, '\t%s (max R2 = %0.4f),\n', n{1}, v(1));
fprintf(1, '\t%s (max R2 = %0.4f), and \n', n{2}, v(2));
fprintf(1, '\t%s (max R2 = %0.4f).\n', n{3}, v(3));

% Best ROI per model
fprintf(1, '\nThe best ROI per predictor is:\n');
for p = 1:numPreds+1
    if p == numPreds+1
        predName = 'Full Model';
        [r2, roi] = max((modelAvgs(:,:,end)), [], 'omitnan');
    else
        predName = predictorList{p};
        [r2, roi] = max((modelAvgs(:,:,end) - modelAvgs(:,:,p)) ./ modelAvgs(:,:,end), [], 'omitnan');
%         [r2, roi] = max(modelAvgs(:,:,end) - modelAvgs(:,:,p), [], 'omitnan');
    end
    roiName = labelStruct(roi).Label;
    fprintf(1,'\t%s: %s (R2 = %0.4f)\n', predName, roiName, r2);
end