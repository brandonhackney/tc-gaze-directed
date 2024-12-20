function val = getMaxMotion(motionTable)
% Given a table of motion vectors over time,
% which are in units of "average pixels displaced per frame",
% identify the maximum value across all elements,
% so we can rescale everything to a more reasonable 0:1
% instead of the current average 1e-4 or whatever.
numStims = height(motionTable);
bigVector = [];
for i = 1:numStims
    bigVector(i) = max(motionTable.MotionEnergy{i});
end
val = max(bigVector);
end