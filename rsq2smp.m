function rsq2smp(data, subNum, varargin)
% Wrapper function to write many SMPs from comparePredictors() output
% Input is a list of model fits per vertex, with two extra dimensions
% data should have i runs, j vertices, and k models tested
% We average across the runs, then generate one SMP per model

sz = size(data);
assert(length(sz) == 3, 'Input 1 must be a 3D matrix!');
numModels = sz(3);

if nargin > 2
    space = varargin{1};
    assert(ischar(space) || isstring(space), 'Third input is a string defining the space (e.g. native, fsaverage)')
else
    space = 'fsaverage';
end

subID = sprintf('sub-%2.f', subNum);
for m = 1:numModels
    label = sprintf('%s_R2-map_model-%i_avg', subID, m);
    writeSMP(100 * mean(data(:, :, m), 1), space, label);
end