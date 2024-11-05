function output = comparePredictors(input)

% First, concatenate all runs for one subject into one big "fullPred"
subNum = 1;
subID = sprintf('sub-%02.f', subNum);
numRuns = 8;
fullPred = [];
dataStack = [];

fprintf(1, 'Getting data for %i runs...', numRuns);
tic;
for r = 1:numRuns
    % Get the design matrix for each run
    tmpPred = getSDM(subNum, r);
    numPredictors = width(tmpPred);
    % May also need to insert nuissance regressors, like head motion
    % Except they're not the same width between runs!!
    nuissance = loadNuissance(subNum, r);
    tmpPred = [tmpPred, nuissance];
    numNuissance = width(nuissance);
    
    % Insert run num as a trailing column
    tmpPred(:, end+1) = r;
    fullPred = [fullPred; tmpPred];
    
    % Also load the actual MRI data for all runs
    dataStack{r} = loadData(subNum, r);
    % Baseline-correct each run independently
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
    else
        % Ignore one of the predictors,
        % to estimate its unique contribution
        pred = fullPred(:,1:width(fullPred) ~= p);
        fprintf(1, 'Predictor %i of %i:\n', p, numPredictors);
    end
    % Iterate through different combinations of data
    iterFits = zeros(numRuns, numVoxels);
    for r = 1:numRuns
        fprintf(1, '\tTesting against run %i/%i...', r, numRuns);
        tic;
        % Do an n-1: caluclate betas for everything BUT r
        % That means we need to subset the data AND the predictors
        trainData = stackRuns(dataStack(1:numRuns ~= r));
        trainPred = pred(pred(:,end) ~= r, :);
        testPred = pred(pred(:,end) == r, :);
        % Deal with run indicators
        trainPred = convertRunCol(trainPred); % conv to many binary cols
        testPred(:,end) = []; % drop the run column for the single run
        % Get the data to be compared
        [betas, residuals] = simpleGLM(trainData, trainPred);
        testData = dataStack{r};
        testResid = averageRunResiduals(residuals, numTRs);
        if p == numPredictors + 1
            % The full model includes one more predictor than "usual"
            % So use the actual numPredictors, not numPredictors - 1
            % This is still different than just betas(:,:),
            % which will include several extra nuissance regressors.
            predictedTS = testPred * betas(1:numPredictors + numNuissance,:);
        else
            predictedTS = testPred * betas(1:numPredictors+numNuissance-1,:);
            % betas above are subset to just the predictors of interest,
            % i.e. it ignores the run number predictors.
        end

        % Get a signed-squared correlation b/w prediction and test data
        % Following McMahon et al 2023
        % This allows for a negative R2, 
        % when the model fits worse than a horizontal line
        iterFits(r,:) = getR2(testData, predictedTS);
%         predCorr(j) = corr2(testData, predictedTS); % j undefined
        % Export to some variable
        fprintf(1, 'Done. ');
        toc % implicitly includes a newline
    end
    
    % After iterating over left-out runs, compare the model fits
    % This... probably depends on the ROI?
    % Just export for now, I'll have to inspect before analyzing further.
%     output(p) = [];
    output{p} = iterFits;
end

% Now after iterating over left-out predictors, compare model fits
% winner = find(max(output));
% if winner == numPredictors + 1
%     fprintf(1, 'Full model is best!');
% else
%     fprintf(1, 'Best to leave out predictor %i', winner);
%     % But... how to find the name of that predictor??
% end