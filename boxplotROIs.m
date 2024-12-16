function boxplotROIs(results, roiLabels, predList)
% Expected input is a 3D matrix of R2s, numSubs * numROIs * numModels
% This visualizes the results as a boxplot, with variance across subs
% Generates a separate figure for every predictor

% Define a list of ROIs we care about:
% V1, MST, FST, MT, 3 STS regions and 3 TPOJ regions
useLabels = [2, 3, 158, 24, 129, 130, 131, 140, 141, 142];

numModels = size(results, 3);
numPreds = length(predList);
if numPreds == numModels
    % Assume the last one is for the full model, and ignore it
    iters = numPreds - 1;
else
    iters = numPreds;
end

% Plot full model first
figure();
    boxplot(results(:,useLabels,end));
    xticklabels({roiLabels(useLabels).Label});
    title('Variance attributable to full model');
    ylabel('R2');
    xlabel('Region of Interest');
    xtickangle(30);
    ylim([0 .06]); % It's unusual to go above 1%, but it may happen.

% Then plot individual predictors
if numPreds > 1
for i = 1:iters
    predName = predList{i};
    figure();
    boxplot(results(:,useLabels,end) - results(:,useLabels,i));
    xticklabels({roiLabels(useLabels).Label});
    title(sprintf('Variance attributable to %s', predName));
    ylabel('R2');
    xlabel('Region of Interest');
    xtickangle(30);
    ylim([0 .06]); % It's unusual to go above 1%, but it may happen.
end
end
    