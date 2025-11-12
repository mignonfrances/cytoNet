% Function to generate a cuttoff correlation value by scrambling each data set, 
% finding all correlation coefficients, and then using the a percentile of the
% set of coefficients as the cutoff value
% INPUT: timeSeries - (nxm) matrix where n is the number of cells and m is
% the number of time points; cutoffPercentile - percentile of scrambled
% correlation coefficients to use as cutoff; M - maximum lag over which to
% calculate cross-correlation values (if M is 0, this becomes the Pearson
% correlation)

function cutoff = generateCutoff(timeSeries, cutoffPercentile, M)
scrambledTimeSeries = timeSeries;
s = size(scrambledTimeSeries);
nCells = s(2);
nFrames = s(1);

% scramble data, while preserving total activity
for j = 1:nCells
    t = randi(nFrames);
    front = scrambledTimeSeries(1:t-1, j);
    back = scrambledTimeSeries(t:nFrames, j);
    index = nFrames - t;
    scrambledTimeSeries(1:(index+1), j) = back;
    scrambledTimeSeries((index+2):nFrames, j) = front;
end

% generate vector of correlation coefficients
% crossCorrelation = [];
% parfor j = 1:nCells
%     scrambledTimeSeries1 = scrambledTimeSeries(:, j);
%     for k = (j+1):nCells
%         scrambledTimeSeries2 = scrambledTimeSeries(:, k);
%         if j ~= k
%             crossCorrelation = [crossCorrelation max(xcov(scrambledTimeSeries1, scrambledTimeSeries2, M, 'coeff'))];
%         end
%     end
% end

A = corr(scrambledTimeSeries);
crossCorrelation = A - diag(diag(A)); % set diagonal elements to zero
crossCorrelation = squareform(crossCorrelation);

% set cutoff
cutoff = percentile(crossCorrelation, cutoffPercentile);
end