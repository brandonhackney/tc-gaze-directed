function [dat] = getStimTiming(subNum, runNum)
% Given the subject NUMBER (not ID), and a run number,
% return data from two files that tells you the stimulus timing

% example fname:
% fname = '/data2/2021_PhysicalSocial/source/SES02/beh/sub-01/sub-01_task-tricopa_date-19-Nov-2021_run-8.txt';
srcID = ind2src(subNum); % use the appropriate subjectID for source
fname = findTxt(srcID, runNum);
dat = readtable(fname, 'Delimiter', '\t');
% The .txt file contains run number, trial number, stim name, 
% and stimulus onset time in both seconds and TRs.
% Does NOT contain stimulus durations.
stimNames = dat.StimID;

% The .tsv and .prt files in the same folder DO contain durations
% So now load one of those up
tsvfname = regexprep(fname, 'date-(\w+)-(\w+)-(\w+)_run', 'run');
tsvfname = strrep(tsvfname, '.txt', '_events.tsv');
tsv = tdfread(tsvfname); % works better than readtable
tsv.stim_id = stimNames; % entire column

% Now stitch the two together
dat.StimName = tsv.stim_id; % StimName was always -1
dat.Duration = tsv.duration;
dat.Offset = tsv.onset + tsv.duration; % tsv.onset is slightly more precise