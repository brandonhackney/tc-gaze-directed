function output = getTrialTiming(subNum)
% Generate an SDM indicating trial timing across all runs for one subject
% Output is a 3D matrix of size n TRs * 2 columns * m trials
% The two columns on each page indicate [this trial, all other trials]
% So column 1 of page 1 indicates when trial 1 happened,
% column 1 of page 2 indicates when trial 2 happened, etc.

numRuns = countRuns(subNum, 'tricopa');
baseTrial = 0; baseTime = 0; baseTR = 0;
bigTrialCol = []; frameCol = []; TRvec = [];
for r = 1:numRuns
    % Load experiment timing data
    dat = getStimTiming(subNum, r);
    % Get timing info at the two sampling rates
    [timeCol, trCol] = trialTimeVectors(subNum, r);
    
    % Get a vector describing trial timing
    numTrials = height(dat);
    numFrames = length(timeCol);
    trialCol = zeros(numFrames, 1);
    for t = 1:numTrials
        onset = dat.Time(t); % sec
        endtime = dat.Offset(t); % sec
        % Those may not line up perfectly with the frame timing, so estimate.
        [~, onsetInd] = min(abs(timeCol - onset));
        [~, offsetInd] = min(abs(timeCol - endtime));
        subset = onsetInd:offsetInd;
        trialCol(subset,1) = t;
    end
    
    % Increment values based on previous runs
    trialCol(trialCol ~= 0) = trialCol(trialCol ~= 0) + baseTrial;
    baseTrial = max(trialCol);
    timeCol = timeCol + baseTime;
    baseTime = max(timeCol);
    trCol = trCol + baseTR;
    baseTR = max(trCol);
    
    % Stack on top of previous runs
    bigTrialCol = [bigTrialCol; trialCol];
    frameCol = [frameCol; timeCol'];
    TRvec = [TRvec; trCol'];
end

% Generate a 3D timing matrix comparing each trial to all others
sdm2 = convertTrialCol(bigTrialCol);
for i = 1:baseTrial
    output(:,:,i) = hrfDownsample(sdm2(:,:,i), frameCol, TRvec);
end

end