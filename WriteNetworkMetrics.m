% Function to write network metrics to file
% INPUT: outputPath - directory path for output, cellInfoAllCells - MATLAB
% structure containing local and global metrics, metric names, adjacency
% matrices, and graph type tag ('spatial' or 'functional') for each image

function [] = writeNetworkMetrics(outputPath, cellInfoAllCells)

globalMetrics = cellInfoAllCells.globalMetrics;
globalMetricNames = cellInfoAllCells.globalMetricNames;
localMetrics = cellInfoAllCells.localMetrics;
localMetricNames = cellInfoAllCells.localMetricNames;
globalMetricsRandom = cellInfoAllCells.globalMetricsRandom;
adjacencyMatrixBinary = cellInfoAllCells.adjacencyMatrixBinary;

nGlobalMetrics = numel(globalMetricNames);
nLocalMetrics = numel(localMetricNames);

fileName = cellInfoAllCells.fileName;
graphTypeTag = cellInfoAllCells.graphTypeTag;

%% print global metrics
writePath = fullfile(outputPath, strcat(fileName, '-GlobalMetrics-', graphTypeTag, '.csv'))
fid = fopen(writePath, 'w');

for i = 1:nGlobalMetrics
    fprintf(fid, '%s,', globalMetricNames{i});
end
fprintf(fid, '\n');

writematrix(globalMetrics, writePath, 'WriteMode', 'append');
fclose(fid);

%% print random graph metrics
writePath = fullfile(outputPath, strcat(fileName, '-GlobalMetricsRandom-', graphTypeTag, '.csv'));
fid = fopen(writePath, 'w');

for i = 1:nGlobalMetrics
    fprintf(fid, '%s,', globalMetricNames{i});
end
fprintf(fid, '\n');

writematrix(globalMetricsRandom, writePath, 'WriteMode', 'append');
fclose(fid);

%% print local metrics
writePath = fullfile(outputPath, strcat(fileName, '-LocalMetrics-', graphTypeTag, '.csv'));
fid = fopen(writePath, 'w');

for i = 1:nLocalMetrics
    fprintf(fid, '%s,', localMetricNames{i});
end
fprintf(fid, '\n');

writematrix(localMetrics, writePath, 'WriteMode', 'append');
fclose(fid);

%% print adjacency matrix
writePath = fullfile(outputPath, strcat(fileName, '-AdjacencyMatrix-', graphTypeTag, '.csv'));
csvwrite(writePath, adjacencyMatrixBinary);

%% save cellInfoAllCells to file
writePath = fullfile(outputPath, strcat(fileName, '-cellInfoAllCells-', graphTypeTag, '.mat'));
save(writePath, 'cellInfoAllCells');
end