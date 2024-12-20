function output = convertTrialCol(trialCol)
% Given a column specifying trial number timing, e.g. 1 1 1 2 2 etc, 
% generate a 3D SDM representing each trial vs all other trials.
% Output is i timepoints by 2 indicators by k trials,
% and trialCol is a vector of i > k timepoints indicating the k trials.

    % Identify unique trials, minus any 0's for non-trial time
    uniqueTrials = unique(trialCol(trialCol ~= 0));
    nTrials = length(uniqueTrials);
    
    % Initialize output matrix
    % Dimensions: i (timepoints) x 2 (binary) x k (trials)
    output = zeros(length(trialCol), 2 , nTrials);
    
    % Loop through each trial
    for t = 1:nTrials
        currentTrial = uniqueTrials(t);
        
        % Logical masks for current trial
        isCurrentTrial = (trialCol == currentTrial);
        isNotCurrentTrial = (trialCol ~= currentTrial & trialCol ~= 0);
        
        % Stack
        output(isCurrentTrial, 1, t) = 1;
        output(isNotCurrentTrial, 2, t) = 1;
    end
end
