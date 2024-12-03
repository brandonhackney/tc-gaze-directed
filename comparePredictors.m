function [output, bics] = comparePredictors(varargin)
% 1. Aggregate MRI data for all runs of a specific task 
% 2. Regress out effects of nuisance parameters (like CSF signal)
% 3. Estimate effect of parameters of interest on the residuals
% 4. Use those betas to predict a timeseries for a left-out run
% 5. Measure the similarity of predicted to actual data
% 6. Cross-validate by leaving out each run independently
% 7. Compare the effect of leaving each predictor out

% First, concatenate all runs for this subject into one big "fullPred"
numRuns = 8;
fullPred = [];
dataStack = [];

hem1 = {'L', 'R'}; % different functions want different formats
hem2 = {'lh', 'rh'};

hem = 2; % 
if nargin > 0
    % Allow a list of subjects to iterate over (even just one)
    subList = varargin{1};
    numSubs = lenth(subList);
else
    numSubs = 30;
    subList = 1:numSubs;
end

% This is hard to preallocate bc I need to know the number of ROIs,
% which is handled by splitByROI() near the end of the pipeline
% size is numSubs * numROIs * numPredictors
% results = zeros([numSubs, 1, 1]); 

for subNum = subList
    % Reset a few things per subject
    fullPred = [];    
    dataStack = [];
    subResults = [];
    % Now go:
    fprintf(1, 'Subject %i: ', subNum);
    fprintf(1, 'Getting data for %i runs...', numRuns);
    tic;
    for r = 1:numRuns
        % Get the design matrix for each run
        [tmpPred, predList, timing] = getSDM(subNum, r);
%         numPredictors = width(tmpPred) + 1;
        numPredictors = length(predList);
        % Also need to account for nuisance regressors, like head motion
        % Except the full files are not the same width between runs!!
        % So only grab a select few predictors
        % And keep these separate, since we will predict in steps
        nuisance = loadNuisance(subNum, r);
%         tmpPred = [tmpPred, nuisance, timing];
        tmpPred = [tmpPred, timing, nuisance];

        % Insert run num as a trailing column
        tmpPred(:, end+1) = r;
        fullPred = [fullPred; tmpPred];

        % Also load the actual MRI data for all runs
        dataStack{r} = loadData(subNum, r, hem1{hem});
        % Baseline-correct each run independently
    %     dataStack{r} = zscore(dataStack{r}, [], 'all');
    %     dataStack{r} = dataStack{r} - mean(dataStack{r}, 1); % subtract mean
    %     dataStack{r} = dataStack{r} ./ max(dataStack{r}, [], 1); % scale to 1
    end
    clear tmpPred
    numVoxels = width(dataStack{1}); % assuming all are same size
    numTRs = height(dataStack{1});
    fprintf(1, 'Done.\n');
    toc

    fprintf(1, '\nCross-validating %i predictors:\n', numPredictors)
    % Now do some leave-one-out analyses
    for p = 1:numPredictors + 1
        % Iterate through different combinations of predictors
        if p == numPredictors+1
            % Use the full model
            pred = fullPred;
            fprintf(1, 'Full model:\n');
            numRunPred = numPredictors;
        else
            % Ignore one of the predictors,
            % to estimate its unique contribution
            pred = fullPred(:,1:width(fullPred) ~= p);
            fprintf(1, 'Predictor %i of %i:\n', p, numPredictors);
            numRunPred = numPredictors - 1;
        end
        % Iterate through different combinations of data
        iterFits = zeros(numRuns, 1); % but will expand...
        iterBICs = zeros(numRuns, 1);
        for r = 1:numRuns
            fprintf(1, '\tTesting against run %i/%i...', r, numRuns);
            tic;
            % Do an n-1: caluclate betas for everything BUT r
            % That means we need to subset the data AND the predictors
            trainData = stackRuns(dataStack(1:numRuns ~= r));
            testData = dataStack{r};
            trainPred = pred(pred(:,end) ~= r, :);
            testPred = pred(pred(:,end) == r, :);
            % Convert run indicators (1 2 3 etc) to an intercept (0/1) per run
            % But since we z-scored the data by run earlier,
            % we don't need to keep these in trainPred:
            % we instead use them to separate nuisance regressors by run
            [trainPred, trainNuis] = splitRunRegressors(trainPred, numRunPred); % conv to many binary cols
    %         testPred(:,end) = []; % drop the run column for the single run
            [testPred, testNuis] = splitRunRegressors(testPred, numRunPred);
            % z-score the nuisance matrix, but not the brain data yet
            trainNuis = zscore(trainNuis);
            testNuis = zscore(testNuis);

            % Also insert an overall intercept column
            trainPred = [ones(height(trainPred), 1), trainPred];
            testPred = [ones(height(testPred), 1), testPred];
            trainNuis = [ones(height(trainNuis), 1), trainNuis];
            testNuis = [ones(height(testNuis), 1), testNuis];
            
            % TRY SUBSETTING RIGHT HERE
            [trainData, trainLabels] = splitByROI(trainData, hem2{hem});
            [testData, testLabels] = splitByROI(testData, hem2{hem});

            % Regress out the nuisance predictors from training data
            [~, trainResid] = simpleGLM(trainData, trainNuis);
            trainResid = zscore(trainResid);
            % Estimate betas for the predictors of interest
            [betas, ~] = simpleGLM(trainResid, trainPred);
            
            % Regress out nuisance from test data
            [~, testResid] = simpleGLM(testData, testNuis);
            testResid = zscore(testResid);
    %         testResid = averageRunResiduals(residuals, numTRs);
            % Calculate expected whole-brain signal based on training model
            predictedTS = testPred * betas; % simple
            
            % Downsample from whole-brain to parcel-average
%             [trainResid, trainLabels] = splitByROI(trainResid, hem2{hem}); % not necessary at this point??
%             [testResid, testLabels] = splitByROI(testResid, hem2{hem});
%             predictedTS = splitByROI(predictedTS, hem2{hem});
            
            % Get a signed-squared correlation b/w prediction and test data
            % Following McMahon et al 2023
            % This allows for a negative R2, 
            % when the model fits worse than a horizontal line
            [x, SSE] = getR2(testResid, predictedTS);
            iterFits(r,1:length(x)) = x; % this helps the variable expand
            % Export to some variable
            iterBICs = BIC(height(testResid), width(testPred), SSE);
            fprintf(1, 'Done. ');
            toc % implicitly includes a newline
        end

        % After iterating over left-out runs, return the model fits
        % This varies by ROI, so we'll need to do an analysis later
        subResults(:,:,p) = iterFits; % where p is the left-out model parameter
        bics(p) = mean(iterBICs);
    end
    % subResults is 8 runs by 181 ROIs by however many parameters
    % average over the runs per subject, then store here:
    results(subNum, :, :) = mean(subResults, 1, 'omitnan');
end % subject

% Now after iterating over left-out predictors, compare model fits
scoreComparison(results, testLabels, predList);

% Generate some QC plots
plotPredCorr(fullPred, subNum, predList);

% EXPORT
% Expand results back from 1-per-parcel to whole-brain
% Necessary for visualization in BrainVoyager
for j = 1:numPredictors + 1
    for i = 1:numRuns        
        output(i,:,j) = expandROIs(results(i,:,j), testLabels);
    end
end
% Write to file
rsq2smp(output, subNum);