function [output, predList, varargout] = getSDM(subNum, runNum, trialStyle)
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


% % WHAT ARE YOU ANALYZING?? % %
predList = {'MotionFrame', 'Interact', 'TopDown', 'Ramp'};
% predList = {'Onset', 'TopDown', 'Interact'};
numPreds = length(predList);
% % WHAT ARE YOU ANALYZING?? % %

% Parse input
if nargin < 3
    trialStyle = false;
end
% BUT
if nargout > 2
    trialStyle = true;
end

% Load data
% example fname:
% fname = '/data2/2021_PhysicalSocial/source/SES02/beh/sub-01/sub-01_task-tricopa_date-19-Nov-2021_run-8.txt';
srcID = ind2src(subNum); % use the appropriate subjectID for source
fname = findTxt(srcID, runNum);
dat = readtable(fname, 'Delimiter', '\t');
% The .txt file contains run number, trial number, stim name, 
% and stimulus onset time in both seconds and TRs.
% Does NOT contain stimulus durations.
numTrials = size(dat, 1);
stimNames = dat.StimID;

% The .tsv and .prt files in the same folder DO contain durations
% So now load one of those up
tsvfname = regexprep(fname, 'date-(\w+)-(\w+)-(\w+)_run', 'run');
tsvfname = strrep(tsvfname, '.txt', '_events.tsv');
tsv = tdfread(tsvfname); % works better than readtable
tsv.stim_id = stimNames; % entire column
% Now you have the onset, duration, and name of each stimulus.

% Get our different predictor tables
params = getPredData(predList);

% The above are at a different sampling rate than the MRI data,
% so we'll need to do quite a bit of math to align them.
% First, specify that the videos have a framerate of 60Hz
SR = 1/60;
% Get the number and duration of TRs for this run:
[numTRs, TR] = findRunLength(subNum, runNum);
durSecs = TR * numTRs;
numFrames = durSecs / SR;
frameCol = 0:SR:durSecs-SR; % use this to look up where to index

% Now generate a predictor matrix:
% rows are timepoints, cols are predictors

sdm = zeros(numFrames, numPreds);
trialCol = zeros(numFrames, 1);
for t = 1:numTrials
    stimName = tsv.stim_id{t};
    if strcmp(stimName(1:2), 'f_')
        % strip out the flip flag and log elsewhere
        % but... what to do with it?
        stimName = stimName(3:end);
        flipFlag = 1;
    else
        flipFlag = 0;
    end
    onset = tsv.onset(t); % sec
    duration = tsv.duration(t); % sec
    endtime = onset + duration; % sec
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
    % Get trial indicator if requested
    if trialStyle
        trialCol(subset,1) = t;
    end
end

% Now you have a boxcar at 60Hz that needs to be:
% - Convolved with an HRF to produce an expected brain response, then
% - Downsampled to the 1.5-sec TR of the MRI data.
% Convolve with the HRF first to preserve as much data as possible,
% relative to downsampling first.
TRvec = 0:TR:(numTRs-1) * TR;
output = hrfDownsample(sdm, frameCol, TRvec);

if trialStyle
    % Generate a 3D timing matrix comparing each trial to all others
    sdm2 = convertTrialCol(trialCol);
    for i = 1:numTrials
        sdm3(:,:,i) = hrfDownsample(sdm2(:,:,i), frameCol, TRvec);
    end
    varargout{1} = sdm3;
end



%% SUBFUNCTIONS
function dataList = buildDataList(predNames, stimName, duration, params)
for p = 1:length(predNames)
    name = predNames{p};
    switch name
        case 'Ramp'
            data = (1:duration) / duration;
        case 'Timing'
            data = ones(duration, 1);
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
        case 'Trial'
            data = ones(duration, 1);
            data = data .* t;
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
end % main function

