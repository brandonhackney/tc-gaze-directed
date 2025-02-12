function [output, predList] = getSDM(subNum, runNum, predList)
% Generate a design matrix for fMRI analysis based on scan and stim data.
% Given a subject number and run number, assuming a specific task,
% find files indicating stimulus order, duration, etc.,
% and generate a statistical design matrix with multiple predictors
% that can be used to calculate betas for that run's brain data.
%
% These folders have .tsv, .mat, and .prt files as well,
% but they are all useles because none indicate stimulus name:
% ONLY the .txt files registered the stimulus name.
% So we need to combine information from two files to get what we need.

global randFlag
if nargin < 3
    % % WHAT ARE YOU ANALYZING?? % %
    predList = {'MotionFrame', 'Interact', 'TopDown', 'Rating'};
%     predList = {'MotionFrame', 'Fixation', 'Interact'};
    % % WHAT ARE YOU ANALYZING?? % %
end
numPreds = length(predList);

% Load experiment timing data
dat = getStimTiming(subNum, runNum);
numTrials = height(dat);
% Now you have the onset, duration, and name of each stimulus.

% Get our different predictor tables
params = getPredData(predList);

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
    % Those may not line up perfectly with the frame timing, so estimate.
%     onsetInd = find(onset <= frameCol, 1);
%     offsetInd = find(endtime >= frameCol, 1, 'last');
    [~, onsetInd] = min(abs(frameCol - onset));
    [~, offsetInd] = min(abs(frameCol - endtime));
    subset = onsetInd:offsetInd;

%     % Get the stuff you need for this stim
%     motion2 = params.motionTable.MotionEnergy{strcmp(motionTable.StimName, stimName)}; % vector
% 
%     % Now start inserting predictors from left to right
%     if length(subset) ~= length(motion2)
%         % Force things to be the same length, in lieu of a better solution
%         fprintf(1, 'Trial %i: subset is %i elements, but motion data is %i elements\n', t, length(subset), length(motion2));
%         subset = onsetInd:(onsetInd + length(motion2) - 1);
%     end
    dataList = buildDataList(predList, stimName, length(subset), params);
    % Insert predictors programmatically
    for s = 1:numPreds
        pname = predList{s};
        x = strcmp({dataList.Name}, pname);
        sdm(subset,s) = dataList(x).Data;
    end
end

if randFlag
    % Shuffle the order of the rows BEFORE applying the HRF
    sdm = sdm(randperm(height(sdm)), :);
end

% Now you have a boxcar at 60Hz that needs to be:
% - Convolved with an HRF to produce an expected brain response, then
% - Downsampled to the 1.5-sec TR of the MRI data.
% Convolve with the HRF first to preserve as much data as possible,
% relative to downsampling first.
output = hrfDownsample(sdm, frameCol, TRvec);

end % main function

%% SUBFUNCTIONS
function dataList = buildDataList(predNames, stimName, duration, params)
for p = 1:length(predNames)
    name = predNames{p};
    switch name
        case 'Ramp'
            data = (1:duration) / duration;
        case 'Timing'
            data = ones(duration, 1);
        case 'Fixation'
            data = params.fixationTable.ScaledFixation(strcmp(params.fixationTable.StimName, stimName)); % scalar
        case 'MotionFrame'
            motion2 = params.motionTable.MotionEnergy{strcmp(params.motionTable.StimName, stimName)}; % vector
            maxMotion = getMaxMotion(params.motionTable);
            data = motion2 ./ maxMotion;
        case 'MotionAvg'
            data = buildDataList({'MotionFrame'}, stimName, duration, params);
            data = mean(data, 1);
        case 'Interact'
            data = params.interactTable2.Interactivity{strcmp(params.interactTable2.StimName, stimName)}; % vector
        case 'InteractAvg'
            data = params.interactTable.Interactivity(strcmp(params.interactTable.StimName, stimName)); % scalar
        case 'Rating'
            data = params.ratingTable.Rating(strcmp(params.ratingTable.StimName, stimName)); % scalar
            maxRating = 5;
            data = data ./ maxRating;
        case 'TopDown'
            data = params.deviatTable.Deviance{strcmp(params.deviatTable.StimName, stimName)}; % vector
            data = interp1(linspace(1,duration,length(data)), data, 1:duration);
            maxDeviation = sqrt(1920^2 + 1200^2);
            data = data ./ maxDeviation;
        case 'TopDownBinary'
            % Instead of the actual deviation values,
            % consider a binary "was this frame deviated?",
            % which could be further reduced to percent time deviated
            data = buildDataList({'TopDown'}, stimName, duration, params);
            data = double(data >= 83.81); % 2 deg visual angle from exp.
            data(data > 1) = 1; % saturation
        case 'Onset'
            data = zeros(duration, 1);
            data(1) = 1;
        case 'Offset'
            data = zeros(duration, 1);
            data(end) = 1;
    end
    
    % Validate size
    if length(data) > duration
        data = data(1:duration);
    elseif length(data) < duration
        d2 = data;
        data = zeros(duration, 1);
        data(1:length(d2)) = d2;
    end
    
    % Write to export
    dataList(p).Name = name;
    dataList(p).Data = data;
end % for predictor
end % subfunction

