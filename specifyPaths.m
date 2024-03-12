function p = specifyPaths(varargin)
% Establishes relative path locations for this project
% An optional input allows you to set a 'base' directory other than pwd

% Establish base directory
if nargin == 0
    p.base = pwd;
else
    assert(isstr(varargin{1}) && isfolder(varargin{1}), 'Input must be a path!');
    p.base = varargin{1};
end

% Define everything else relative to p.base
% Assuming directory structure is BIDS/derivatives/here:
p.derivatives = fullfile(p.base, '..');
p.bv = fullfile(p.derivatives, 'bv'); % BrainVoyager-format data
p.fs = fullfile(p.derivatives, 'sourcedata', 'freesurfer');

end