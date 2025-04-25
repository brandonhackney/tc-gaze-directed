% Generate a t-Map using the AQ-discriminative videos as a contrast
% i.e. if AQ discriminates this video, it's ToM-heavy. If not, ToM-light.
% 

subNum = 1; numRuns = 8;
fullPred = [];
for r = 1:numRuns
    % Get the design matrix for each run
    [tmpPred, predList] = getSDMContrast(subNum, r);
    tmpPred = zscore(tmpPred);
    % Insert an intercept column
    tmpPred = [ones(height(tmpPred), 1), tmpPred];
    % Insert run num as a trailing column, for later subsetting
    tmpPred(:, end+1) = r;
    fullPred = [fullPred; tmpPred]; % stack runs
end


% Extract data from gifti files and stack all runs into one big matrix
% [dataStack, roiLabels] = getDataStack(subNum, 2);
dataStack = getDataStack(subNum, 2, 'fsnative');
data = stackRuns(dataStack);
contrast = [0, 1, -1, 0]; % constant, c1, c2, runNum. % convert run to cols
[~,~,~,tMap] = simpleGLM(data, fullPred, 0, contrast);
tMap = tMap(:); % ensure it's vertical

% Write data to a gifti
out = gifti();
if numel(tMap) < 100000
    % Use roiLabels to write the tMap to surface
    cdata = zeros(size(roiLabels(1).Verts));
    for i = 1:length(roiLabels)
        v = roiLabels(i).Verts;
        cdata(v) = tMap(i);
    end
    out.cdata = cdata;
else
    out.cdata = tMap;
end

% Threshold
% hack temp value is just look at t beyond +/-3
out.cdata(out.cdata < 3 & out.cdata > -3) = 0;

% Visualize on surface
fname2 = '/data2/2021_PhysicalSocial/BIDS/derivatives/sub-01/ses-01/anat/sub-01_ses-01_hemi-R_inflated.surf.gii';
surf = gifti(fname2);
figure(); plot(surf, out); colorbar;
title(sprintf('%s ToM-heavy > ToM-light', subID));