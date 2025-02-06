function [frameCol, TRvec] = trialTimeVectors(subNum, runNum)
% First, specify that the videos have a framerate of 60Hz
SR = 1/60;
% Get the number and duration of TRs for this run:
[numTRs, TR] = findRunLength(subNum, runNum);
% Calculate the timing vectors
durSecs = TR * numTRs;
frameCol = SR:SR:durSecs; % use this to look up where to index
TRvec = TR:TR:(numTRs * TR); % use this to get the time of each TR
% use the END time of the TR instead of the start time,
% because the start time would ignore anything that happens DURING the TR