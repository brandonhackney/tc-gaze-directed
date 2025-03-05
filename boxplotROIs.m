function boxplotROIs(results, roiLabels, predList)
% Expected input is a 3D matrix of R2s, numSubs * numROIs * numModels
% This visualizes the results as a boxplot, with variance across subs
% Generates a separate figure for every predictor

% Define a list of ROIs we care about:
% V1, V4, V4t, V6, MST, FST, MT, 3 STS regions, 3 TPOJ regions, FEF, LO1-3
useLabels = [2, 7, 157, 4, 3, 158, 24, 129, 130, 131, 140, 141, 142, 11, 21, 22, 160];
useLabels = ismember(1:length(roiLabels), useLabels); % convert to logical
% useLabels = 1:length(roiLabels); % use all
% useLabels = max(results(:,:,end),[],1) > .0484; % if any sub above threshold
% useLabels = mean(results(:,:,end), 1, 'omitnan') > .0484; % threshold picked via permutation testing
numUsed = sum(useLabels);

numModels = size(results, 3);
numPreds = length(predList);
if numPreds < numModels
    % Assume the last one is for the full model, and ignore it
    iters = numModels - 1;
    rtype = 'R2';
else
    iters = numPreds;
    rtype = 'r';
end

% Preprocess roi names
for i = 1:length(roiLabels)
    x = roiLabels(i).Label;
    x = erase(x, '_ROI');
    x = x(3:end); % erase 'R_xxx' part
    roiLabels(i).Label = x;
end

if strcmp(rtype, 'R2')
% Plot full model first
figure();
    plot(0:numUsed+1, zeros(size(0:numUsed+1)), '--'); % draw a line at 0
    hold on;
    boxplot(results(:,useLabels,end));
    hold off;
    xticklabels({roiLabels(useLabels).Label});
    title('Variance attributable to full model');
    ylabel('R2');
    xlabel('Region of Interest');
    xtickangle(30);
    ylim([0 1]); % 0-100%. It's unusual to go above 1%, but it may happen.
end

% Then plot individual predictors
if numPreds > 1
for i = 1:iters
    predName = predList{i};
    figure();
    if strcmp(rtype, 'R2')
        plot(0:numUsed+1, zeros(size(0:numUsed+1)), '--'); % draw a line at 0
        hold on;
        boxplot((results(:,useLabels,end) - results(:,useLabels,i)));
        hold off;
        title(sprintf('Unique variance attributable to %s', predName));
        ylabel('R2_f - R2_i');
        ylim([-.05 .15]);
%         ytickformat('percentage');
    elseif strcmp(rtype, 'r')
        plot(0:numUsed+1, zeros(size(0:numUsed+1)), '--'); % draw a line at 0
        hold on;
        boxplot(results(:,useLabels, i));
        hold off;
        title(sprintf('Correlation of %s with timeseries', predName));
        ylabel('Pearson''s r');
        ylim([-1 1]);
    end
    xticklabels({roiLabels(useLabels).Label});
    xlabel('Region of Interest');
    xtickangle(30);
    
end

% Make another figure grouping by ROI, across predictors
warning('off', 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout');
figure();
tiledlayout(ceil(sqrt(numUsed) / 2), ceil(sqrt(numUsed) * 2));
for i = find(useLabels)
    nexttile;
    plot(0:numModels + 1, zeros(size(0:numModels + 1)), '--');
    hold on;
    title(strrep(roiLabels(i).Label, '_', '\_'));
    if strcmp(rtype, 'R2')
        f0 = squeeze(results(:,i,end));
        f1 = squeeze(results(:,i,1:end-1));
        boxplot(f0 - f1);
        ylabel('R2_f - R2_i');
        ylim([-.05 .15]);
        yticks(-.05:.025:.15);
    else
        boxplot(results(:,i,:));
        ylabel('Pearson''s r with timeseries');
        ylim([-1 1]);
    end
    hold off;
    xticklabels(predList);
    xtickangle(30);
end
warning('on', 'MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout');
end