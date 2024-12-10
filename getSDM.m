function [output, predList, timing] = getSDM(subNum, runNum)
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
predList = {'MotionFrame', 'Interact', 'TopDown'};
numPreds = length(predList);
% % WHAT ARE YOU ANALYZING?? % %

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
motionTable = importdata('motionDataSum.mat'); % total, not avg
% motionTable = importdata('motionData.mat');
interactTable = importdata('interactData.mat');
ratingTable = importdata('rateData.mat');
interactTable2 = importdata('interactVectors.mat');
deviatTable = importdata('devData.mat');
% ...

% The above are at a different sampling rate than the MRI data,
% so we'll need to do quite a bit of math to align them.
% First, specify that the videos have a framerate of 60Hz
SR = 1/60;
TR = 1.5; % but try to read this in from a file somewhere
% Get the number of TRs for this run:
numTRs = findRunLength(subNum, runNum);
durSecs = TR * numTRs;
numFrames = durSecs / SR;
frameCol = 0:SR:durSecs-SR; % use this to look up where to index

% Now generate a predictor matrix:
% rows are timepoints, cols are predictors

sdm = zeros(numFrames, numPreds);
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
    onsetInd = find(onset <= frameCol, 1);
    offsetInd = find(endtime >= frameCol, 1, 'last');
    subset = onsetInd:offsetInd;

    TR = dat.Time(2) / dat.TR(2); % units of TRs (expect 1.5 sec)
    timeVec = onset:TR:endtime;
    TRvec = timeVec / TR; % should basically be 1 2 3 etc
    col = ones([length(timeVec),1]); % a prototype  

    % Get the stuff you need for this stim
    motion2 = motionTable.MotionEnergy{strcmp(motionTable.StimName, stimName)}; % vector

    % Now start inserting predictors from left to right
    if length(subset) ~= length(motion2)
        % Force things to be the same length, in lieu of a better solution
        subset = onsetInd:(onsetInd + length(motion2) - 1);
    end
    dataList = buildDataList(predList, length(subset));
    % Insert predictors programmatically
    for s = 1:numPreds
        pname = predList{s};
        x = strcmp({dataList.Name}, pname);
        sdm(subset,s) = dataList(x).Data;
    end
end

% Experimental: address collinearity of parametric predictors
% Stimulus timing is highly correlated with all other predictors,
% especially when you just modulate its amplitude for different trials.
% so mean-center (NOT z-score!) anything that isn't timing.
x = ~strcmp(predList, 'Timing') & ~strcmp(predList, 'Ramp');
sdm(:,x) = sdm(:,x) - mean(sdm(:,x));
sdm(:,x) = sdm(:,x) ./ max(sdm(:,x));

% Now you have a boxcar at 60Hz that needs to be:
% - Convolved with an HRF to produce an expected brain response
% - Downsampled to the 1.5-sec TR of the MRI data.
% Convolve with the HRF first to preserve as much data as possible,
% relative to downsampling first.
hrf = spm_hrf(SR);

TRvec = 0:TR:(numTRs-1) * TR; % 
for i = 1:width(sdm)
    col = sdm(:,i);
    col = conv(col, hrf, 'same');
    % Vectors need to be summed/averaged over the longer time period
%     mtvec = (1:1000*SR:length(col) * 1000*SR) - 1;
%     output(:,i) = binData(frameCol', col, TRvec'); % average the framewise values within each TR
    output(:,i) = interp1(frameCol, col, TRvec);
end

if nargout > 2 && any(contains(predList, 'Timing'))
    % Return the timing column (trial on/off) as a separate variable.
    tind = strcmp(predList, 'Timing');
    timing = output(:,tind);
    output(:,tind) = [];
elseif nargout > 2 && ~any(contains(predList, 'Timing'))
    timing = []; % give it an empty to avoid crashing
end

% SUBFUNCTIONS
function dataList = buildDataList(predNames, duration)
for p = 1:length(predNames)
    name = predNames{p};
    switch name
        case 'Ramp'
            data = (1:duration) / duration;
        case 'Timing'
            data = 1;
        case 'MotionFrame'
            motion2 = motionTable.MotionEnergy{strcmp(motionTable.StimName, stimName)}; % vector
            maxMotion = getMaxMotion(motionTable);
            data = motion2 ./ maxMotion;
        case 'MotionAvg'
            data = buildDataList({'MotionFrame'}, duration);
            data = mean(data, 1);
        case 'Interact'
            data = interactTable2.Interactivity{strcmp(interactTable2.StimName, stimName)}; % vector
        case 'InteractAvg'
            data = interactTable.Interactivity(strcmp(interactTable.StimName, stimName)); % scalar
        case 'Rating'
            data = ratingTable.Rating(strcmp(ratingTable.StimName, stimName)); % scalar
            maxRating = 5;
            data = data ./ maxRating;
        case 'TopDown'
            data = deviatTable.Deviance{strcmp(deviatTable.StimName, stimName)}; % vector
            data = interp1(linspace(1,duration,length(data)), data, 1:duration);
            maxDeviation = sqrt(1920^2 + 1200^2);
            data = data ./ maxDeviation;
        case 'TopDownBinary'
            % Instead of the actual deviation values,
            % consider a binary "was this frame deviated?",
            % which could be further reduced to percent time deviated
            data = buildDataList({'TopDown'}, duration);
            data = double(data >= 83.81); % 2 deg visual angle from exp.
            data(data > 1) = 1; % saturation
    end
    % Write to export
    dataList(p).Name = name;
    dataList(p).Data = data;
end % for predictor
end % subfunction
end % main function

