function labelStruct = scoreComparison(modelFits, labelStruct, predictorList)
% Run from within comparePredictors
% Take the smaller stack of model fits (i.e. one per ROI, not whole-brain)
% Return the best model (and associated R2) per ROI
numModels = size(modelFits, 3);
numPreds = length(predictorList);
numROIs = size(modelFits, 2);

% Get the average model fits across subjects,
% so we can see which model has the best fit in each ROI.
modelAvgs = mean(modelFits, 1);

% Now you have a 1xn vector of model numbers per ROI
% Write that, and the associated R2 (or average R2 I guess), to the struct
assert(length(labelStruct) == numROIs, 'Model fits and labels have unequal number of ROIs!');
roimeans = zeros(numROIs, 1);
roimaxes = zeros(numROIs, 1);
for i = 1:numROIs
    results = zeros(1,numPreds);
    % Variance partitioning:
    % R2 of full model MINUS R2 of model missing one predictor
    % Reports how much each predictor contributes to the overall regression
    fprintf(1, '%s:\t',labelStruct(i).Label);
    fprintf(1, 'Variance attributable to:');
    for p = 1:numPreds+1
        if p == numPreds + 1
            fprintf('\t%s: %0.4f', 'FullModel', modelAvgs(:,i,end));
%             fprintf('\t%s: %0.4f', 'FullModel', 1); % meaningless if x-y/x
        else
            predName = predictorList{p};
            % Do the variance partitioning within subject, then average
            fullmodel = modelFits(:,i,end);
            limmodel = modelFits(:,i,p);
            result = (fullmodel - limmodel) ./ fullmodel;
            results(p) = mean(result, 1);
            fprintf('\t%s: %0.4f', predName, results(p));
        end
    end
    fprintf(1, '\n');
    
    % Write these results into a struct too
    [r2, x] = max(results);
    labelStruct(i).BestModel = predictorList{x};
    labelStruct(i).R2 = r2; % also export R2 val
    
    % Measure consistency of each ROI
    roimeans(i) = mean(modelFits(:,i,end));
    roimaxes(i) = max(modelFits(:,i,end));
end

fprintf(1, '\n\nMean Max = %0.4f\tMax Mean = %0.4f\n', mean(roimaxes, 'omitnan'), max(roimeans));
fprintf(1, 'Max Max = %0.4f\tMean Mean = %0.4f\n', max(roimaxes), mean(roimeans, 'omitnan'));

% Best-explained ROIs (i.e. peak average across all subjects)
[v, i] = maxk(roimeans, 3);
n = {labelStruct(i).Label};
fprintf(1, '\nThe three most consistently well-explained ROIs are:\n');
fprintf(1, '\t%s (mean R\x00B2 = %0.4f),\n', n{1}, v(1));
fprintf(1, '\t%s (mean R\x00B2 = %0.4f), and \n', n{2}, v(2));
fprintf(1, '\t%s (mean R\x00B2 = %0.4f).\n', n{3}, v(3));

% Single-subject best model fit
[v, i] = maxk(roimaxes, 3);
n = {labelStruct(i).Label};
fprintf(1, '\nThe three highest single-subject fits are in:\n');
fprintf(1, '\t%s (max R\x00B2 = %0.4f),\n', n{1}, v(1));
fprintf(1, '\t%s (max R\x00B2 = %0.4f), and \n', n{2}, v(2));
fprintf(1, '\t%s (max R\x00B2 = %0.4f).\n', n{3}, v(3));

% Best ROI per model
fprintf(1, '\nThe best ROI per predictor is:\n');
for p = 1:numPreds+1
    if p == numPreds+1
        predName = 'Full Model';
        r2 = 1;
        [r2f, roi] = max((modelAvgs(:,:,end)), [], 'omitnan');
    else
        predName = predictorList{p};
        fullmodel = modelFits(:,:,end);
        limmodel = modelFits(:,:,p);
        result = (fullmodel - limmodel) ./ fullmodel;
        result = mean(result, 1, 'omitnan'); % collapse subjects
        
        [r2, roi] = max(result, [], 'omitnan');
%         [r2, roi] = max(modelAvgs(:,:,end) - modelAvgs(:,:,p), [], 'omitnan');
        r2f = modelAvgs(:,roi,end);
    end
    roiName = labelStruct(roi).Label;
    fprintf(1,'\t%s: %s (%0.2f%% of R\x00B2 = %0.4f)\n', predName, roiName, r2 * 100, r2f);
end