function data = loadData(sub, run, varargin)

% Parse inputs
if nargin > 2
    hemi = varargin{1};
    assert(ischar(hemi) || isstring(hemi), 'Third input must be a string indicating hemisphere!');
    assert(length(hemi) <=2, 'Hemisphere should be simply ''L'' or ''R''');
    hemi = char(upper(hemi(1))); % force it to conform to expectations
else
    hemi = 'R';
end

% Generate a filename and find the full path to it
pths = specifyPaths;
subID = sprintf('sub-%02.f', sub);
name = ['hemi-', hemi, '_space-fsaverage_bold.func.gii'];
fname = findSubData(pths, subID, run, name);

% Load gifti data using NIH's plugin
assert(exist('gifti', 'file') == 2, 'Cannot find gifti plugin! Download from https://github.com/gllmflndn/gifti')
data = gifti(fname);
% Result is an object with a single field: cdata
% data.cdata is an n*t matrix of n vertices by t timepoints
% If you care about the locations of the vertices at all,
% that data comes from... elsewhere.
% For now, extract and rotate 90 deg to meet expectations of existing functions:
data = data.cdata'; % overwrite to save memory
end