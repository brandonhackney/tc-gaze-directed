%% Setup
% Currently just loading the first run
% Will eventually want to concatenate multiple runs
pths = specifyPaths;
sub = 'sub-01';
run = 1;
name = 'hemi-R_space-fsaverage_bold.func.gii';
fname = findSubData(pths, sub, run, name);

%% Pattern
% Load gifti data using NIH's plugin
assert(exist('gifti', 'file') == 2, 'Cannot find gifti plugin! Download from https://github.com/gllmflndn/gifti')
data = gifti(fname);
% Result is an object with a single field: cdata
% data.cdata is an n*t matrix of n vertices by t timepoints
% If you care about the locations of the vertices at all,
% that data comes from... elsewhere.
% For now, extract and rotate 90 deg to meet expectations of existing functions:
data = data.cdata'; % overwrite to save memory

%% Predictors
% To calculate betas, you need a predictor matrix and a contrast.
% Pred should be n timepoints * p predictors, including constant
predfname = ['*.sdm']; % need an actual name
% See null_batchGLM line 43:56 for code to handle multiple runs
file = xff(predfname);
pred = file.SDMMatrix;
file.clearobj;
%% Contrast
% This was copied verbatim from the atlas project.
% It may be inappropriate for the current analysis, 
% since we're not really doing contrasts per se.
conditionList = importdata('getFilePartsFromContrast.mat');
contName = conditionList(t).contrast;
[posInd,negInd, ~] = getConditionFromFilename(contName);
contrast = getContrastVector(size(pred,2),posInd,negInd);

%% Calculate t-statistic map
% This function can also generate whole-brain betas and residuals,
% but we'll skip that for now.
[tMap,~,~] = simpleGLM(data,pred,contrast);