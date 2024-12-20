function output = getPredData(predList)
% Variably load predictor data depending on which ones you're using.
% This really only gets called by one function, but hey.

    output = struct();
    if any(strcmp(predList, 'MotionFrame')) || any(strcmp(predList, 'MotionAvg'))
        output.motionTable = importdata('motionDataSum.mat'); % total, not avg
    end
    % motionTable = importdata('motionData.mat');
    if any(strcmp(predList, 'InteractAvg'))
        output.interactTable = importdata('interactData.mat');
    end
    if any(strcmp(predList, 'Rating'))
        output.ratingTable = importdata('rateData.mat');
    end
    if any(strcmp(predList, 'Interact'))
        output.interactTable2 = importdata('interactVectors.mat');
    end
    if any(contains(predList, 'TopDown'))
        output.deviatTable = importdata('devData.mat');
    end
    % ...add any new metrics here
    
end % function