function [output, bics] = comparePredictors(subNum)
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

fprintf(1, 'Getting data for %i runs...', numRuns);
tic;
for r = 1:numRuns
    % Get the design matrix for each run
    tmpPred = getSDM(subNum, r);
    numPredictors = width(tmpPred);
    % Also need to account for nuisance regressors, like head motion
    % Except the full files are not the same width between runs!!
    % So only grab a select few predictors
    % And keep these separate, since we will predict in steps
    nuisance = loadNuisance(subNum, r);
    tmpPred = [tmpPred, nuisance];
    
    % Insert run num as a trailing column
    tmpPred(:, end+1) = r;
    fullPred = [fullPred; tmpPred];
    
    % Also load the actual MRI data for all runs
    dataStack{r} = loadData(subNum, r, hem1{hem});
    % Baseline-correct each run independently
    dataStack{r} = zscore(dataStack{r}, [], 'all');
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
        trainPred = pred(pred(:,end) ~= r, :);
        testPred = pred(pred(:,end) == r, :);
        % Convert run indicators (1 2 3 etc) to an intercept (0/1) per run
        % But since we z-scored the data by run earlier,
        % we don't need to keep these in trainPred:
        % we instead use them to separate nuisance regressors by run
        [trainPred, trainNuis] = splitRunRegressors(trainPred, numRunPred); % conv to many binary cols
%         testPred(:,end) = []; % drop the run column for the single run
        [testPred, testNuis] = splitRunRegressors(testPred, numRunPred);
        % Also insert an overall intercept column
        trainPred = [ones(height(trainPred), 1), trainPred];
        testPred = [ones(height(testPred), 1), testPred];
        % Regress out the nuisance predictors from training data
        [~, trainResid] = simpleGLM(trainData, trainNuis);
        
        % Subset the data to parcels (from whole-brain)
        [trainResid, trainLabels] = splitByROI(trainResid, hem2{hem});
        
        % Estimate betas for the predictors of interest
        [betas, ~] = simpleGLM(trainResid, trainPred);
        % Regress out nuisance from test data
        testData = dataStack{r};
        [~, testResid] = simpleGLM(testData, testNuis);
        % Subset the test data to parcels 
        [testResid, testLabels] = splitByROI(testResid, hem2{hem});
%         testResid = averageRunResiduals(residuals, numTRs);
        predictedTS = testPred * betas; % simple

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
    
    % After iterating over left-out runs, compare the model fits
    % This... probably depends on the ROI?
    % Just export for now, I'll have to inspect before analyzing further.
%     output(p) = [];
    results(:,:,p) = iterFits;
    bics(p) = mean(iterBICs);
end

% Now after iterating over left-out predictors, compare model fits
% winner = find(max(output));
% if winner == numPredictors + 1
%     fprintf(1, 'Full model is best!');
% else
%     fprintf(1, 'Best to leave out predictor %i', winner);
%     % But... how to find the name of that predictor??
% end

% EXPORT
% Expand results back from 1-per-parcel to whole-brain
% Necessary for visualization in BrainVoyager
for j = 1:numPredictors + 1
    for i = 1:numRuns        
        output(i,:,j) = expandROIs(results(i,:,j), testLabels);
    end
end