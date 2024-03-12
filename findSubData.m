function output = findSubData(pths, subID, run, input)
% Locate an individual file inside this mess of folders
% Requires the output of specifyPaths as first input
% Then give it a subject ID (e.g. 'sub-03'),
% a run number (e.g. [2]), 
% and the rest of your desired filename
% (e.g. 'hemi-R_space-fsaverage_bold.func.gii')
% I should probably do more about that last part...

% Validate input
assert(isstruct(pths), 'First input was not a struct! Expected output of specifyPaths()');
assert(ischar(subID), 'Second input was not a string! Expected subject ID');
assert(isnumeric(run), 'Third input was non-numeric! Expected a run number');
assert(ischar(input), 'Fourth input was not a string! Expected a partial filename');

% Locate data
datdir = pths.derivatives;
fpath = fullfile(datdir, subID, 'ses-02', 'func');
S = {subID, 'ses-02', 'task-tricopa', ['run-', run], input};
fname = strjoin(S, '_');
output = fullfile(fpath, fname);

% Validate output before exiting
assert(exist(output, 'file'), 'Desired file %s not found', output);
end