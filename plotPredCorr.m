function plotPredCorr(predCorr, predLabels)
% Plot a correlation matrix as a heatmap, with variable labels.

% Define the expected list of predictors. (probably should be an input).
% predLabels = {'Motion', 'Interactivity', 'Top-Down Gaze', 'Understanding', 'StimTiming'};
numPred = length(predLabels);
assert(width(predCorr) >= numPred, 'Not enough columns in this predictor matrix!');

predLabels = strrep(predLabels, '_', '\_');

figure();
h = heatmap(predCorr);
h.XDisplayLabels = predLabels;
h.YDisplayLabels = predLabels;

if numPred > 10
    axp = struct(h);
    set(axp.NodeChildren(3), 'XTickLabelRotation', 60);
    set(axp.NodeChildren(3), 'YTickLabelRotation', 0);
    h.FontSize = 6.5;
    axp.Axes.XAxisLocation = 'top';
end

title('Predictor collinearity');
xlabel('Averaged across all subjects');

% Scale colormap for correlations (bounded -1 to +1)
if exist('clim', 'file')
    clim([-1, 1]);
else
    caxis([-1, 1]);
end
colormap('turbo');

end