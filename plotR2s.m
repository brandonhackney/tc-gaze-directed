function plotR2s(r2Stack, subNum)
% Generate a subplot window with histograms showing whole-brain R2s per run
% Intended to run for a single subject at a time: second input is sub num
% Assume data orientation is rows == run, cols == vertices

numRuns = height(r2Stack);
figure();
tl = tiledlayout(1,numRuns);
for i = 1:numRuns
    nexttile;
    histogram(r2Stack(i,:));
    title(sprintf('Run %i', i));
    xlim([0 1]);
end
title(tl, sprintf('Subject %i', subNum));