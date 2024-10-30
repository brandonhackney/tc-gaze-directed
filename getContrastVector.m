function contvec = getContrastVector(numPred,posInd,negInd)
% contvec = getContrastVector(numPred,posInd,negInd)
% numPred is the number of predictors in a design matrix (ie num of betas)
% posInd is a vector of indices for the positive conditions
% negInd is a vector of indices for the negative conditions
% contvec is a scaled vector defining the contrast based on the inputs
% e.g. if numPred = 4, posInd = [1 2], and negInd = 3, contvec = [1 1 -2 0]

contvec = zeros([1,numPred]);
contvec(posInd) = length(negInd);
contvec(negInd) = -1 * length(posInd);
end