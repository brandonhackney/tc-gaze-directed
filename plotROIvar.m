function plotROIvar(modelFits, labelStruct, predList)
% After running comparePredictors, send me the R2s and the ROI labels
% I will plot the variance attributable to each predictor

i = 2; % 24 for MT, 4 for V6

roiName = labelStruct(i).Label;
roiFits = modelFits(:,i,:); % row = subject, col = roi, page = model
roiFits = squeeze(roiFits); % row = subject, col = model

ttxt = sprintf('%s\nError bars are over subjects', roiName);
ttxt = strrep(ttxt, '_', '\_');

figure();
boxplot(roiFits * 100);
xticklabels([predList, {'Full Model'}]);
xtickangle(30);
ylabel('ROI variance explained');
ytickformat('percentage');
title(ttxt);


roiProp = (roiFits - roiFits(:,end)) ./ roiFits(:,end);
figure();
boxplot(roiProp * 100);
xticklabels([predList, {'Full Model'}]);
xtickangle(30);
ylabel('Proportion of full-model variance explained');
title(ttxt);
ylim([-200, 200]);
ytickformat('percentage');
end