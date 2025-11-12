% Script calling functions to generate functional graph, generate
% functional network plots, compute graph-based metrics, and write metrics
% to file
% INPUT: errorReport = struct containing information on errors during processing, 
% filePath - directory path to input file, assumed to be a stacked
% .tif file containing calcium imaging data; outputDirName - directory path
% to folder for output; optionalMask - optional binary mask file to define
% ROIs

% DEPENDENCIES:
% Mike Wu (2021). wgPlot - Weighted Graph Plot (a better version of gplot)
% https://www.mathworks.com/matlabcentral/fileexchange/24035-wgplot-weighted-graph-plot-a-better-version-of-gplot
% MATLAB Central File Exchange

function errorReport = calciumEngine(errorReport, filePath, outputDirName, maskPath)

[~, fileName, ~] = fileparts(filePath);


% generate functional graph based on cross-correlation analysis
if nargin < 4
    [errorReport, cellInfoAllCells] = generateFunctionalGraph(errorReport, filePath);
else
    [errorReport, cellInfoAllCells] = generateFunctionalGraph(errorReport, filePath, maskPath);
end

if sum(cellInfoAllCells.activeCells) == 0
    errorReport(end+1) = makeErrorStruct(['no active cells found in file ', fileName], 1);
    return;
end

% compute graph-based metrics
[cellInfoAllCells.globalMetrics, cellInfoAllCells.localMetrics, ...
    cellInfoAllCells.globalMetricNames, cellInfoAllCells.localMetricNames, processParams] = ...
    calculateNetworkMetrics(cellInfoAllCells.adjacencyMatrixBinary);
writeLog('[calciumEngine] network metrics calculated');

% compute network metrics for random graph
cellInfoAllCells.globalMetricsRandom = calculateNetworkMetricsRandom(cellInfoAllCells);
writeLog('[calciumEngine] network metrics for random graph calculated');

if ~processParams.processCompleted
    errorReport(end+1) = makeErrorStruct(['processing for file ', fileName, ' exceeds computing capacity; graph not created'], 1);
end

% generate processed image plots and write metrics to file
generateFunctionalNetworkPlots(cellInfoAllCells, fileName, outputDirName);
writeLog('[calciumEngine] processed plots generated');
writeNetworkMetrics(outputDirName, cellInfoAllCells);
writeLog('[calciumEngine] metrics written to file');
end