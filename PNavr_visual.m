% Select multiple recording folders
recordingDirs = {};
defaultPath = pwd;  % Start from current folder

while true
    thisDir = uigetdir(defaultPath, 'Select a recording folder (Cancel to finish)');
    if thisDir == 0
        break;  % Cancel pressed
    end

    % Check if all 41 step folders contain 'clamp1.ma'
    numSteps = 41;
    allPresent = true;

    for stepIdx = 0:numSteps-1
        stepStr = sprintf('%03d', stepIdx);
        fileToCheck = fullfile(thisDir, stepStr, 'Clamp1.ma');
        if ~isfile(fileToCheck)
            allPresent = false;
            break;
        end
    end
    
  defaultPath = thisDir;  % Update default path for convenience
  
    if allPresent
        recordingDirs{end+1} = thisDir;
       
    else
        % Show warning dialog
        msg = sprintf('Recording folder "%s" is incomplete.\n\nMissing Clamp1.ma file in one or more step folders (000–040).\n\nIt will be excluded.', thisDir);
        uiwait(warndlg(msg, 'Incomplete Folder Skipped', 'modal'));
    end
end


% Number of steps (assumed fixed)
numSteps = 41;

% Preallocate cell array to hold data from each step
stepData = cell(numSteps, 1);

% Loop through each step
for stepIdx = 0:numSteps-1
    stepStr = sprintf('%03d', stepIdx);
    stepTraces = {};

    % Loop through each selected recording
    for r = 1:length(recordingDirs)
        stepPath = fullfile(recordingDirs{r}, stepStr, 'Clamp1.ma');

        % Load the data (replace with your preferred loading method if needed)
        if isfile(stepPath)
            trace = loadClamp1MA(stepPath);  % See helper function below
            stepTraces{end+1} = trace;
        else
            warning('Missing: %s', stepPath);
        end
    end

    % Average across recordings for this step
    if ~isempty(stepTraces)
        traceMat = cell2mat(cellfun(@(x) x(:), stepTraces, 'UniformOutput', false));
        avgTrace = mean(traceMat, 2);

        % Store averaged trace
        stepData{stepIdx+1} = avgTrace;
    end
end

% Ask where to save averaged result
saveDir = uigetdir(pwd, 'Select folder to save averaged .ma steps');

% Create a figure to hold all overlaid averaged traces
overlayFig = figure('Name', 'Overlaid Averaged Traces', 'NumberTitle', 'off');
hold on;
colors = lines(numSteps);  % Distinct color for each trace


% Save each averaged step
for stepIdx = 0:numSteps-1
    stepStr = sprintf('%03d', stepIdx);
    stepFolder = fullfile(saveDir, stepStr);
    if ~exist(stepFolder, 'dir')
        mkdir(stepFolder);
    end

   avgTrace = stepData{stepIdx+1};
if ~isempty(avgTrace)
    % Save averaged trace as .ma file
    data = avgTrace;
    save(fullfile(stepFolder, 'Clamp1.ma'), 'data');

    % Add trace to overlay figure
    plot(avgTrace, 'LineWidth', 1.2, 'Color', colors(stepIdx+1, :), 'DisplayName', ['Step ', stepStr]);
end
end

    % Finalize and save the overlay figure
xlabel('Sample');
ylabel('Voltage (mV)');
title('Overlaid Averaged Traces (All Steps)');
legend('show');
grid on;
hold off;

% Save the overlay figure
overlayFigPath = fullfile(saveDir, 'AllSteps_Overlay.png');
saveas(overlayFig, overlayFigPath);

% Show the figure interactively (optional)
uiwait(overlayFig);  % Wait for user to close

disp('Averaging complete and .ma files saved.');

%% Helper function
function data = loadClamp1MA(filePath)
    try
        allData = h5read(filePath, '/data');
        
        % Check shape and select only voltage trace
        if size(allData, 1) == 3
            data = allData(1, :);  % Only the voltage trace
        elseif size(allData, 2) == 3
            data = allData(:, 1);  % Also voltage trace
        elseif size(allData, 3) == 3
            data = allData (:, 1);  % Also voltage trace
        else
            data = allData;  % Assume already just voltage
        end
    catch ME
        error('Failed to load %s: %s', filePath, ME.message);
    end
end
