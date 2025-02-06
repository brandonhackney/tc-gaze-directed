function labelStruct = scoreComparison(modelFits, labelStruct, predictorList)
% Run from within comparePredictors
% Take the smaller stack of model fits (i.e. one per ROI, not whole-brain)
% Return the best model (and associated R2) per ROI
[~, numROIs, numModels] = size(modelFits);
numPreds = length(predictorList);

if numModels > numPreds
    predictorList{end+1} = 'FullModel';
    numPreds = numPreds + 1;
    rtype = 'R2';
else
    rtype = 'r';
end

% Get the average model fits across subjects,
% so we can see which model has the best fit in each ROI.
modelAvgs = mean(modelFits, 1);

% Now you have a 1xn vector of model numbers per ROI
% Write that, and the associated R2 (or average R2 I guess), to the struct
assert(length(labelStruct) == numROIs, 'Model fits and labels have unequal number of ROIs!');
roimeans = zeros(numROIs, 1);
roimaxes = zeros(numROIs, 1);
results = zeros(numROIs, numPreds);
for i = 1:numROIs
    roiResults = zeros(1,numPreds);
    % Variance partitioning:
    % R2 of full model MINUS R2 of model missing one predictor
    % Reports how much each predictor contributes to the overall regression
    fprintf(1, '%s:\t',labelStruct(i).Label);
    fprintf(1, 'Variance attributable to:');
    for p = 1:numPreds
        predName = predictorList{p};
        limmodel = modelFits(:,i,p);
        if strcmp(rtype, 'R2')
            % R-squares imply that we want to do variance partitioning
            % Do it within subject, then average
            if strcmp(predName, 'FullModel')
                % Don't compare the full model to itself;
                % just report the overall R2
                result = modelAvgs(:,i,end);
            else
                % Compare each predictor to the full model
                fullmodel = modelFits(:,i,end);
                result = (fullmodel - limmodel) ./ fullmodel;
            end
        elseif strcmp(rtype, 'r')
            % These are correlations, so just report them straight off
            % "limited model" is a misnomer in this case.
            result = limmodel;
        end
        roiResults(p) = mean(result, 1);
        fprintf('\t%s: %0.4f', predName, roiResults(p));
    end
    results(i,:) = roiResults;
    fprintf(1, '\n');
    if strcmp(rtype, 'R2')
        % Drop the full-model results
        roiResults(end) = [];
    end
    % Write these results into a struct too
    [r2, x] = max(roiResults);
    labelStruct(i).BestModel = predictorList{x};
    labelStruct(i).R2 = r2; % also export R2 val
    
    % Measure consistency of each ROI
    roimeans(i) = mean(modelFits(:,i,end));
    roimaxes(i) = max(modelFits(:,i,end));
end

if strcmp(rtype, 'R2')
    % This block analyzes the full models.
    % There are no equivalent metrics for correlations, so skip.   
    fprintf(1, '\n\nMean Max = %0.4f\tMax Mean = %0.4f\n', mean(roimaxes, 'omitnan'), max(roimeans));
    fprintf(1, 'Max Max = %0.4f\tMean Mean = %0.4f\n', max(roimaxes), mean(roimeans, 'omitnan'));
    
    % Best-explained ROIs (i.e. peak average across all subjects)
    [v, i] = maxk(roimeans, 3);
    n = {labelStruct(i).Label};
    fprintf(1, '\nThe three most consistently well-explained ROIs are:\n');
    fprintf(1, '\t%s (mean R\x00B2 = %0.4f),\n', n{1}, v(1));
    fprintf(1, '\t%s (mean R\x00B2 = %0.4f), and \n', n{2}, v(2));
    fprintf(1, '\t%s (mean R\x00B2 = %0.4f).\n', n{3}, v(3));
end

% Single-subject best model fit
[v, i] = maxk(roimaxes, 3);
n = {labelStruct(i).Label};
if strcmp(rtype, 'R2')
    fprintf(1, '\nThe three highest single-subject fits are in:\n');
    fprintf(1, '\t%s (max R\x00B2 = %0.4f),\n', n{1}, v(1));
    fprintf(1, '\t%s (max R\x00B2 = %0.4f), and \n', n{2}, v(2));
    fprintf(1, '\t%s (max R\x00B2 = %0.4f).\n', n{3}, v(3));
elseif strcmp(rtype, 'r')
    fprintf(1, '\nThe three highest single-subject fits are in:\n');
    fprintf(1, '\t%s (max r = %0.4f),\n', n{1}, v(1));
    fprintf(1, '\t%s (max r = %0.4f), and \n', n{2}, v(2));
    fprintf(1, '\t%s (max r = %0.4f).\n', n{3}, v(3));
end

% Best ROI per model
fprintf(1, '\nThe best ROI per predictor is:\n');
for p = 1:numPreds
    predName = predictorList{p};
    if strcmp(predName, 'FullModel')
        r2 = 1;
        [r2f, roi] = max((modelAvgs(:,:,end)), [], 'omitnan');
    else
%         [r2, roi] = max(abs(results(:,p)), [], 'omitnan');
        % Which ROI is best explained by this predictor?
        % We don't just want the highest unique variance in isolation:
        % Explaining 100% of 0.0000001% is meaningless.
        % Let's find the highest overall variance explained:
        % the product of the full-model R2 and the predictor's percentage.
        % 5% of 100% is better than 100% of 1%.
        % ...but also, 100% of 0.9% is maybe better than 1% of 100%??
        [~, roi] = max(results(:,p) .* modelAvgs(:,:,end)', [], 'omitnan');
        r2 = results(roi, p);
        r2f = modelAvgs(:,roi,end);
    end
    roiName = labelStruct(roi).Label;
    if strcmp(rtype, 'R2')
        fprintf(1,'\t%s: %s (%0.2f%% of R\x00B2 = %0.4f)\n', predName, roiName, r2 * 100, r2f);
    elseif strcmp(rtype, 'r')
        fprintf(1,'\t%s: %s (Average r = %0.4f)\n', predName, roiName, r2);
    end
end