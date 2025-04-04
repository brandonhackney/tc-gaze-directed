function makeNewPRT(subNum, runNum)
% Generate a PRT to contrast AQ-sensitive stims with the others.
% Relies on a list of stims where eyetracking/rating correlates with AQ

% Load experiment timing data
dat = getStimTiming(subNum, runNum);
% Now you have the onset, duration, and name of each stimulus.

% Now load a list of videos defining condition 1 (implies condition 2)
% These are videos where clarity ratings are most modulated by AQ score
load('sigVids.mat', 'sigVids');

% Compare and determine which vids in dat belong to which condition
cond1 = ismember(dat.StimID, sigVids);
cond2 = ~cond1;

% Open output file
pths = specifyPaths;
fdir = pths.bv;
fname = sprintf('sub-%02.f_task-tricopa_style-contrast_run-%i.prt', subNum, runNum);
fpath = fullfile(fdir, fname);
fid = fopen(fpath, 'w+');

% Write data to file
% Header
fprintf(fid, ['\n',...
'FileVersion:\t\t2\n',...
'\n',...
'ResolutionOfTime:\tmsec\n',...
'\n',...
'Experiment:\t\ttricopa\n',...
'\n',...
'BackgroundColor:\t0 0 0\n',...
'TextColor:\t\t255 255 255\n',...
'TimeCourseColor:\t255 255 255\n',...
'TimeCourseThick:\t3\n',...
'ReferenceFuncColor:\t0 0 80\n',...
'ReferenceFuncThick:\t2\n',...
'\n',...
'NrOfConditions:\t\t2\n',...
'\n',...
]);

% Condition 1: ToM-hard videos
fprintf(fid, 'ToM-hard\n%i\n', sum(cond1));
fprintf(fid, '%8.0f\t%8.0f\n', [round(dat.Time(cond1)' .* 1000); round(dat.Duration(cond1)' .* 1000)]);
fprintf(fid, 'Color:\t\t\t 255 0 0\n');
fprintf(fid, '\n');

% Condition 2: ToM-easy videos
fprintf(fid, 'ToM-easy\n%i\n', sum(cond2));
fprintf(fid, '%8.0f\t%8.0f\n', [round(dat.Time(cond2)' .* 1000); round(dat.Duration(cond2)' .* 1000)]);
fprintf(fid, 'Color:\t\t\t0 255 0\n');
fprintf(fid, '\n');

% Close file
fclose(fid);
