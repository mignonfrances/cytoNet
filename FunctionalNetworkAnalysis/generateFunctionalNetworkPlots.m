function errorReport = generateFunctionalNetworkPlots(cellInfoAllCells, fileName, outputDirName)

errorReport = [];
fileNameBase = createBaseName(fileName);

%% correlation between distance and functional edge strength
cellLocationsConnectedCells = cellInfoAllCells.cellLocations(find(cellInfoAllCells.activeCells), :);
cellDistancesConnectedCells = squareform(pdist(cellLocationsConnectedCells, 'Euclidean'));
avgeDistanceConnectedCells = mean(cellDistancesConnectedCells(cellInfoAllCells.adjacencyMatrixWeighted~=0));

f1 = figure('Visible', 'Off');
set(gcf, 'Color', 'w');
if ~isempty(cellDistancesConnectedCells)
    plot(cellDistancesConnectedCells, cellInfoAllCells.adjacencyMatrixWeighted, 'k.', 'MarkerSize', 14);
end

if isnan(cellInfoAllCells.cutoffCorrelation)
    cellInfoAllCells.cutoffCorrelation = 0;
end

title(sprintf('Average Distance between connected cells = %4.3f pixels', avgeDistanceConnectedCells));
if ~isnan(cellInfoAllCells.cutoffCorrelation)
    r = refline(0, cellInfoAllCells.cutoffCorrelation);
    r.Color = 'r';
end

xlabel('Distance (pixels)');
ylabel('Functional Edge Strength');
set(gca, 'FontSize', 14);

% Save
if ~exist(outputDirName, 'dir')
    mkdir(outputDirName);
end

invalidChars = ['\', '/', ':', '*', '?', '"', '<', '>', '|'];
fileNameBase = regexprep(fileNameBase, ['[', regexptranslate('escape', invalidChars), ']'], '_');
savePath = fullfile(outputDirName, strcat(fileNameBase, '-AveCorrDistance.svg'));

disp(['Saving average correlation distance figure to: ', savePath]);

set(f1, 'Visible', 'on');
drawnow;  % ensure it renders
saveas(f1, savePath);

close(f1);

%% plot weighted graph on maxImage (only correlations above cutoff)
f2 = figure('Visible', 'Off');
set(gcf, 'Color', 'w');
imshow(cellInfoAllCells.maxImage); hold on;
maskBoundaries = bwboundaries(cellInfoAllCells.mask);
for k = 1:length(maskBoundaries)
    currentBoundary = maskBoundaries{k};
    plot(currentBoundary(:, 2), currentBoundary(:, 1), 'r');
end

for k = 1:cellInfoAllCells.nNodes
    text(cellInfoAllCells.cellLocations(k, 1), cellInfoAllCells.cellLocations(k, 2), num2str(k), 'color', 'r', 'FontSize', 8);
end

A = cellInfoAllCells.adjacencyMatrixWeighted;
A(cellInfoAllCells.adjacencyMatrixWeighted < cellInfoAllCells.cutoffCorrelation) = 0;
wgPlot(A, cellInfoAllCells.cellLocations(find(cellInfoAllCells.activeCells), :), ...
    'edgeColorMap', parula);
set(findall(gcf,'type','line'), 'LineWidth', 1);

% Save
if ~exist(outputDirName, 'dir')
    mkdir(outputDirName);
end

invalidChars = ['\', '/', ':', '*', '?', '"', '<', '>', '|'];
savePath = fullfile(outputDirName, strcat(fileNameBase, '-DistanceCorrelation.svg'));

disp(['Saving distance correlation figure to: ', savePath]);

set(f2, 'Visible', 'on');
drawnow;
saveas(f2, savePath);
close(f2);

%% heatmap of number of spikes overlaid on mask

f3 = plotMetricOnImage(cellInfoAllCells.cellLocations, cellInfoAllCells.mask, cellInfoAllCells.nSpikes);
hold on;

for k = 1:cellInfoAllCells.nNodes
    text(cellInfoAllCells.cellLocations(k, 1), cellInfoAllCells.cellLocations(k, 2), num2str(k), 'color', 'k', 'FontSize', 8);
end

savePath = fullfile(outputDirName, strcat(fileNameBase, 'NumberSpikes.svg'));
disp(['Saving spike count heatmap figure to: ', savePath]);
drawnow;
saveas(f3, savePath);
set(f3,'Visible','on');
close(f3);

end