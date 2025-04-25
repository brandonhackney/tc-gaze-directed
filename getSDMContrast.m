function [output, predList] = getSDMContrast(subNum, runNum)

load('sigVids.mat', 'sigVids');
numPreds = 2; % on/off
predList = {'ToM-heavy', 'ToM-light'};

% Load experiment timing data
dat = getStimTiming(subNum, runNum);
numTrials = height(dat);
% Now you have the onset, duration, and name of each stimulus.
% The parameters are at a different sampling rate than the MRI data,
% so we'll need to do quite a bit of math to align them.
% First, get some timing vectors at the different sampling rates
[frameCol, TRvec] = trialTimeVectors(subNum, runNum);
numFrames = length(frameCol);

% Now generate a predictor matrix:
% rows are timepoints, cols are predictors
sdm = zeros(numFrames, numPreds);
for t = 1:numTrials
    stimName = dat.StimName{t};
    if strcmp(stimName(1:2), 'f_')
        % strip out the flip flag and log elsewhere
        % but... what to do with it?
        stimName = stimName(3:end);
        flipFlag = 1;
    else
        flipFlag = 0;
    end
    
    onset = dat.Time(t); % sec
    duration = dat.Duration(t); % sec
    endtime = dat.Offset(t); % sec

    [~, onsetInd] = min(abs(frameCol - onset));
    [~, offsetInd] = min(abs(frameCol - endtime));
    subset = onsetInd:offsetInd;
    
    % Condition logic:
    if ismember(stimName, sigVids)
        sdm(subset,1) = 1;
    else
        sdm(subset,2) = 1;
    end

end

% Now you have a boxcar at 60Hz that needs to be:
% - Convolved with an HRF to produce an expected brain response, then
% - Downsampled to the 1.5-sec TR of the MRI data.
% Convolve with the HRF first to preserve as much data as possible,
% relative to downsampling first.
output = hrfDownsample(sdm, frameCol, TRvec);

end % main function
