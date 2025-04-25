function [output, roiLabels, predList] = comparePredictors(hem, varargin)
% 1. Aggregate MRI data for all runs of a specific task 
% 2. Regress out effects of nuisance parameters (like CSF signal)
% 3. Estimate effect of parameters of interest on the residuals
% 4. Use those betas to predict a timeseries for a left-out run
% 5. Measure the similarity of predicted to actual data
% 6. Cross-validate by leaving out each run independently
% 7. Compare the effect of leaving each predictor out
global randFlag
predList = {'MotionFrame', 'Interact', 'TopDown', 'Rating'};

hem1 = {'L', 'R'}; % different functions want different formats
hem2 = {'lh', 'rh'};

if nargin < 1
    % Default hem is RH
    hem = 2; 
else
    assert(isnumeric(hem), 'Input 1 must be a number 1 or 2 (Left or Right hemisphere)');
    assert(hem < 3, 'Input 1 must be a number 1 or 2 (Left or Right hemisphere)');
end
numSubs = 0;
if nargin > 1
    % Allow a list of subjects to iterate over (even just one)
    subList = varargin{1};
    numSubs = length(subList);
end
if numSubs == 0
    % if no input OR if empty input
    numSubs = 30;
    subList = 1:numSubs;
end
if nargin > 2
    randFlag = true;
    numIter = varargin{2};
    assert(isnumeric(numIter), 'Second input must be a number of random iterations to run');
    assert(floor(numIter) == numIter, 'Second input must be an integer number of random iterations to run');
    rng('shuffle');
else
    randFlag = false;
    numIter = 1;
end

% This is hard to preallocate bc I need to know the number of ROIs,
% which is handled by splitByROI() near the end of the pipeline
% size is numSubs * numROIs * numPredictors
% results = zeros([numSubs, 1, 1]); 

predCorr = [];
for subNum = subList
    % First, concatenate all runs for this subject into one big stack
    [dataStack, roiLabels] = getDataStack(subNum, hem);
    % Same for predictors
    [fullPred, predList] = getPredictorStack(subNum, predList);
    numPredictors = length(predList);
    numRuns = height(dataStack);
    
%     % Chop out any rows that were dropped from the predictor matrix
%     % i.e. stimuli we are ignoring
%     for r = 1:numRuns
%         dataStack{r}(chop{r},:) = 0;
%     end
    
    % Analyze predictor collinearity
    predCorr(:,:,subNum) = corr(fullPred(:, 2:end-1), 'rows', 'complete');
    disp(predList);
    disp(predCorr(:,:,subNum));
    
    for i = 1:numIter
        if randFlag
            itime = tic;
            % Shuffle the rows of the predictor matrix for each iteration
            fullPred = getPredictorStack(subNum, predList);
        end
        fprintf(1, '\nCross-validating %i predictors:\n', numPredictors)
        % Now do some leave-one-out analyses
        for p = 1:numPredictors + 1
            % Iterate through different combinations of predictors
            if p == numPredictors+1
                % Use the full model
                pred = fullPred;
                fprintf(1, 'Full model:\n');
                numRunPred = numPredictors;
                usedPreds = predList;
            else
                % Ignore one of the predictors,
                % to estimate its unique contribution
                % Use p+1 to account for the intercept column in col 1
                spec = 1:width(fullPred) ~= p+1;
                pred = fullPred(:,spec);
                fprintf(1, 'Predictor %i of %i:\n', p, numPredictors);
                numRunPred = numPredictors - 1;
                usedPreds = predList(spec(2:end-1)); % skip intercept & run
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

                % Drop the run indicator from the predictors.
                % We've already z-scored the data per run,
                % so additional controls are unnecessary.
                testPred(:,end) = [];
                trainPred(:,end) = [];

                % Estimate betas for the predictors of interest
                [betas, ~] = simpleGLM(trainData, trainPred, 1);

                % Calculate expected whole-brain signal based on training model
                predictedTS = testPred * betas; % simple

                % Measure model fit b/w prediction and test data
                [R2, SSE] = getR2(testData, predictedTS);

                iterFits(r,1:length(R2)) = R2; % this helps the variable expand
                % Export to some variable
                iterBICs = BIC(height(testData), width(testPred), SSE);
                fprintf(1, 'Done. ');
                toc % implicitly includes a newline
            end

            % After iterating over left-out runs, return the model fits
            % This varies by ROI, so we'll need to do an analysis later
            % subResults is (run, ROI, model)
            subResults(:,:,p) = iterFits; % where p is the left-out model parameter
            bics(p) = mean(iterBICs);
        end
        if randFlag
            % Collapse the results for this iteration... somehow.
            % I guess just take the overall mean? But preserve predictors
%             thresh(subNum, i,:) = mean(subResults, [1,2], 'omitnan');
            % NO: export the full-model R2, collapsing run, keeping ROI
            % thresh is then subNum, iteration, ROI
            thresh(subNum, i, :) = mean(subResults(:,:,end), 1, 'omitnan');
            fprintf(1, 'Done with iteration %i of %i. ', i, numIter);
            toc(itime);
        end
        
        
    end % for random iterations
    % subResults is 8 runs by 181 ROIs by however many parameters
    % average over the runs per subject, then store here:
    results(subNum, :, :) = mean(subResults, 1, 'omitnan');
end % subject

if randFlag
    try
        % Determine a selection threshold based on the 95th percentile of r2s.
        % thresh is the r2 for the full model only, org as (sub iter ROI)
        % we want the simplest threshold possible,
        % so instead of exporting an expected r2 for each ROI,
        % we'll average across the whole brain. But not as a first step.
        % Average over subjects first,
        output = squeeze(mean(thresh, 1, 'omitnan'));
        % then sort each ROI independently and take the value at position 95%,
        % considering it's now (iter, ROI), default sort() behavior works.
        output = sort(output);
        critInd = floor(.95 * numIter);
            if critInd < 1
                % unlikely scenario but
                critInd = 1;
                fprintf(1, 'Warning! Only 1 iteration performed??\n')
            end
        % then average across ROIs to get your final whole-brain threshold.
        output = mean(output(critInd, :), 'omitnan');
        fprintf(1, '\n\nExpected R2 =\t%0.4f\n', output);
    catch ME
        output = thresh;
        ME
    end
else


% Now after iterating over left-out predictors, compare model fits
scoreComparison(results, roiLabels, predList);
boxplotROIs(results,roiLabels,predList);

% Generate some QC plots
avgPredCorr = mean(predCorr, 3);
plotPredCorr(avgPredCorr, predList);

% EXPORT
% Expand results back from 1-per-parcel to whole-brain
% Necessary for visualization in BrainVoyager
% for j = 1:numPredictors + 1
%     for i = 1:numRuns        
%         output(i,:,j) = expandROIs(results(i,:,j), testLabels);
%     end
% end
% % Write to file
% rsq2smp(output, subNum);
output = results;
end