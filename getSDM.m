function [output, timing] = getSDM(subNum, runNum)
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
numPreds = 4 + 1; % 4 of interest, plus timing as a nuisance variable
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
    motion = mean(motion2, 1); % scalar
    interact = interactTable.Interactivity(strcmp(interactTable.StimName, stimName)); % scalar
    interact2 = interactTable2.Interactivity{strcmp(interactTable2.StimName, stimName)}; % vector
    rating = ratingTable.Rating(strcmp(ratingTable.StimName, stimName)); % scalar
    deviation2 = deviatTable.Deviance{strcmp(deviatTable.StimName, stimName)}; % vector
    
    % Associated parameters
    maxRating = 5;
    maxDeviation = 1920 * 1200;
    maxMotion = getMaxMotion(motionTable);

    % Now start inserting predictors from left to right
    if length(subset) ~= length(motion2)
        % Force things to be the same length, in lieu of a better solution
        subset = onsetInd:(onsetInd + length(motion2) - 1);
    end
%     sdm(subset,1) = 1; % constant
%     sdm(subset,2) = t; % trial
%     sdm(subset,3) = motion;
%     sdm(subset,4) = interact;
%     sdm(subset,5) = rating / maxRating; % convert to percent

%     sdm(subset, 1) = motion; % avg motion across whole video
    sdm(subset,1) = motion2 ./ maxMotion;
    sdm(subset,2) = interact2; % Binary - no need to rescale
    sdm(subset,3) = deviation2 ./ maxDeviation;
    sdm(subset,4) = rating / maxRating; % convert to percent
    sdm(subset,5) = 1; % stim on vs stim off
    % ...
%     sdm = [sdm; dmat];
end
% sdm(:,1) = 1; % constant

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

% Pull the final column (trial on/off) out and return as a separate var.
timing = output(:,end);
output(:,end) = [];
end % function