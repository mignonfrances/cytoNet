pctRunOnAll addpath(genpath(fullfile(matlabroot, 'toolbox', 'stats')));

% Ensure correct pdist function is in use
p = which('pdist', '-all');
if contains(p{1}, ['toolbox' filesep 'stats' filesep 'eml'])
    warning('Removing internal EML path for pdist compatibility.');
    rmpath(fileparts(p{1}));
    rehash toolboxcache;
end

errorReport = struct('errorMessage', {}, 'errorSeverity', {});
filePath = "path_to_your_.tif_file";
outputDirName = "path_to_folder_where_outputs_will_be_saved";
%maskPath = 'path_to_mask';  % Optional, depending on the function's requirements

% Print time 
fprintf('Trial start time %s\n', datestr(now,'HH:MM:SS.FFF'));

% Call the function
errorReport = calciumEngine(errorReport, filePath, outputDirName);
% errorReport = calciumEngine(errorReport, filePath, outputDirName, maskPath) % if using a mask 

% Check the error report
disp(errorReport);

% Print time 
fprintf('Trial finished run at time %s\n', datestr(now,'HH:MM:SS.FFF'));