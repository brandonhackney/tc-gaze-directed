pths = specifyPaths;
sub = 'sub-01';
run = 1;
name = 'hemi-R_space-fsaverage_bold.func.gii';
fname = findSubData(pths, sub, run, name);

% Load gifti data using NIH's plugin
assert(exist('gifti', 'file') == 2, 'Cannot find gifti plugin! Download from https://github.com/gllmflndn/gifti')
data = gifti(fname);
% Result is an object with a single field: cdata
% data.cdata is an n*t matrix of n vertices by t timepoints
% If you care about the locations of the vertices at all,
% that data comes from... elsewhere.
