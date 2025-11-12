% Function to extract calcium time series data from an image stack and
% compute a functional graph based on correlations among time series
%
% INPUT: filePath - stacked tif image file you wish to generate a graph for;
% mask (optional) - binary mask defining ROIs for time series extraction;
% if no mask is provided, the default segmentation program is run on
% the maximum intensity projection image
%
% OUTPUT: cellInfoAllCells - MATLAB structure containing calcium image
% analysis parameters

function [errorReport, cellInfoAllCells] = generateFunctionalGraph(errorReport, filePath, maskPath)

%% Setting up
nSpikesThreshold = 2; % minimum number of spikes to determine an active cell
cutoffPercentile = 99; % percentile of scrambled correlations to be used for cutoff
maxLag = 0; %  maximum value of lag allowed for calculating correlations

% fileName = getFileName(filePath);
[~, fileName, ~] = fileparts(filePath);
infoImage = imfinfo(filePath);

maxImage = findMaxImage(filePath); % extract maximum intensity projection image

% if mask not provided as input, run default segmentation
if nargin < 3
    mask = segmentWatershed(maxImage);
else
    infoMask = imfinfo(maskPath);
    mask = imread(maskPath);
    % reshape mask to fit size of calcium image stack
    if infoMask.Height ~= infoImage(1).Height || infoMask.Width ~= infoImage(1).Width
        mask = imresize(mask, [infoImage(1).Height, infoImage(1).Width]);
    end
end

totalFrames = size(infoImage, 1);
CC = bwconncomp(mask);
nCells = CC.NumObjects;

if nCells == 0
    errorReport(end+1) = makeErrorStruct(['no cells found in file ', fileName], 1);
    cellInfoAllCells = [];
    return;
end

%% gathering information on ROIs
cellInfoAllCells = struct;
cellInfoAllCells.mask = mask;
cellInfoAllCells.maxImage = maxImage;
cellInfoAllCells.adjacencyType = 3; % setting adjacencyType as 3 for functional graphs
cellInfoAllCells.graphTypeTag = 'functional';
cellInfoAllCells.fileName = fileName;

stats = regionprops(CC, 'Centroid');

% centroid locations of all cells
cellInfoAllCells.cellLocations = zeros(nCells, 2);
for k = 1:nCells
    current = stats(k);
    location = current.Centroid;
    cellInfoAllCells.cellLocations(k, 1) = location(1); % storing x location
    cellInfoAllCells.cellLocations(k, 2) = location(2); % storing y location
end

% store unprocessed time-varying intensity for each ROI
cellInfoAllCells.Fraw = zeros(totalFrames, nCells);
for j = 1:totalFrames
    currentFrame = imread(filePath, j);
    stats = regionprops(CC, currentFrame, 'PixelValues', 'Area');
    intensitiesCurrentFrame = zeros(1, nCells);
    parfor k = 1:nCells
        current = stats(k);
        intensitiesCurrentFrame(k) = sum(current.PixelValues)/current.Area;
    end
    cellInfoAllCells.Fraw(j, :) = intensitiesCurrentFrame;
end

writeLog('[generateFunctionalGraph] data extracted');

%% Process data: calculate F-F0/F0

cellInfoAllCells.Fprocessed = zeros(totalFrames, nCells);
cellInfoAllCells.F0 = zeros(totalFrames, nCells);
Fraw = cellInfoAllCells.Fraw;

    for j = 1:nCells % loop over all cells
        currentF = Fraw(:, j);

        % subtract linear trend from mean intensity data
        currentF = detrend(currentF) + currentF(1);

        % create moving baseline vector
        window = round(totalFrames/100); % moving window is a hundredth of total time
        currentF = smooth(currentF, window); % smooth time-varying intensity trace
        F0 = zeros(size(currentF));
        parfor k = 1:totalFrames
            if k <= window
                F0(k) = percentile(currentF(1:k+window), 8);
            elseif k >= length(currentF) - window
                F0(k) = percentile(currentF(k-window:length(currentF)), 8);
            else
                F0(k) = percentile(currentF(k-window:k+window), 8);
            end
        end
    
    cellInfoAllCells.F0(:, j) = F0; % store baseline intensity
    cellInfoAllCells.Fprocessed(:, j) = smooth((currentF - F0)./F0); % calculate delta F / F, smooth data, and store as Fprocessed
end

writeLog('[generateFunctionalGraph] data processed');

%% Filter out inactive cells to avoid false positives and measure activity
cellInfoAllCells.activeCells = zeros(1, nCells); % marker vector for active cells
cellInfoAllCells.nSpikes = zeros(1, nCells); % storing number of spikes for each cell
for j = 1:nCells
    currentSignal = cellInfoAllCells.Fprocessed(:, j);
    
    % A spike is a peak whose deltaF/F value is at least 0.25 of the max value in the signal, or 0.05, whichever is higher
    [spikes, ~] = findpeaks(currentSignal, 'MinPeakHeight', max(0.25*max(currentSignal), 0.05), 'MinPeakProminence', max(0.25*max(currentSignal), 0.05));
    
    cellInfoAllCells.nSpikes(j) = numel(spikes); % storing number of spikes
    
    % setting marker vector to 1 if number of spikes exceeds threshold
    if numel(spikes) > nSpikesThreshold
        cellInfoAllCells.activeCells(j) = 1;
    end
end

%% Generate adjacency matrix based on cross-correlation of pixel intensities, only for active cells
nNodes = sum(cellInfoAllCells.activeCells);
timeSeries = cellInfoAllCells.Fprocessed(:, find(cellInfoAllCells.activeCells));

cutoff = generateCutoff(timeSeries, cutoffPercentile, maxLag); % calculate cutoff using scrambled data
cellInfoAllCells.cutoffCorrelation = cutoff;
writeLog('[generateFunctionalGraph] cutoff generated');

A = corr(timeSeries);
adjacencyMatrixWeighted = A - diag(diag(A)); % set diagonal elements to zero

% create binary adjacency matrix
adjacencyMatrixBinary = ones(size(adjacencyMatrixWeighted));
adjacencyMatrixBinary(adjacencyMatrixWeighted < cutoff) = 0;
cellInfoAllCells.adjacencyMatrixBinary = adjacencyMatrixBinary;
cellInfoAllCells.adjacencyMatrixWeighted = adjacencyMatrixWeighted;
cellInfoAllCells.nNodes = nNodes;
writeLog('[generateFunctionalGraph] adjacency matrix created');
end