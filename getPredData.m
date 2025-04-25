function output = getPredData(predList)
% Variably load predictor data depending on which ones you're using.
% This really only gets called by one function, but hey.
% Using load() is a lot faster than importdata() when you know what's there

    output = struct();
    if any(strcmp(predList, 'MotionFrame')) || any(strcmp(predList, 'MotionAvg'))
        load('motionDataSum.mat', 'motion'); % total, not avg
        output.motionTable = motion;
    end
    % motionTable = importdata('motionData.mat');
    if any(strcmp(predList, 'InteractAvg'))
        load('interactData.mat', 'intScore');
        output.interactTable = intScore;
    end
    if any(strcmp(predList, 'Rating'))
        load('rateData.mat', 'ratings');
        output.ratingTable = ratings;
    end
    if any(strcmp(predList, 'Interact'))
        load('interactVectors.mat', 'intVectors');
        output.interactTable2 = intVectors;
    end
    if any(contains(predList, 'TopDown'))
        load('devData.mat', 'deviance');
        output.deviatTable = deviance;
    end
    if any(contains(predList, 'Fixation'))
        load('fixationData.mat', 'fixations');
        output.fixationTable = fixations;
    end
    if any(contains(predList, 'ToM'))
        load('tomData.mat', 'ToM');
        output.tomTable = ToM;
        output.tomTable.CommBeta = zscore(output.tomTable.CommBeta);
    end
    % ...add any new metrics here
    
end % function