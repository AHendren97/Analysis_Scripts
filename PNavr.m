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

% Preallocate cell array: one cell per step, each containing valid traces
stepData = cell(numSteps, 1);
stepTracesAll = cell(numSteps, 1);  % Each cell will hold traces from accepted recordings

% Preallocate cell array: one cell per step, each containing valid traces
stepData = cell(numSteps, 1);
stepTracesAll = cell(numSteps, 1);  % Each cell will hold traces from accepted recordings

% Loop through each recording (ask user once per recording)
for r = 1:length(recordingDirs)
    recordingPath = recordingDirs{r};

    % Try to load all steps
    thisRecTraces = cell(numSteps, 1);
    valid = true;

    for stepIdx = 0:numSteps-1
        stepStr = sprintf('%03d', stepIdx);
        stepPath = fullfile(recordingPath, stepStr, 'Clamp1.ma');

        if isfile(stepPath)
            trace = loadClamp1MA(stepPath);
            thisRecTraces{stepIdx+1} = trace;
        else
            warning('Missing step %s in %s', stepStr, recordingPath);
            valid = false;
            break;
        end
    end

    % Skip this recording if incomplete
    if ~valid
        continue;
    end

    % Plot all traces together
    hFig = figure('Name', ['Recording ', num2str(r), ' (All Steps)'], 'NumberTitle', 'off');
    hold on;
    for stepIdx = 0:numSteps-1
        plot(thisRecTraces{stepIdx+1}, 'DisplayName', sprintf('Step %03d', stepIdx));
    end
    xlabel('Sample'); ylabel('Voltage (mV)');
    title(['Recording ', num2str(r), ' - All Steps']);
    legend('show'); grid on;

    % Ask once: include this recording?
    choice = questdlg(sprintf('Include recording %d in averaging?', r), ...
        'Recording Review', ...
        'Yes', 'No', 'Yes');

    close(hFig);  % Close plot window

    if strcmp(choice, 'Yes')
        % Append each trace to step-wise pool
        for stepIdx = 0:numSteps-1
            stepTracesAll{stepIdx+1}{end+1} = thisRecTraces{stepIdx+1};
        end
    else
        fprintf('Excluded: %s\n', recordingPath);
    end
end

% Now average across included recordings for each step
for stepIdx = 0:numSteps-1
    traces = stepTracesAll{stepIdx+1};
    if ~isempty(traces)
        traceMat = cell2mat(cellfun(@(x) x(:), traces, 'UniformOutput', false));
        avgTrace = mean(traceMat, 2);
        stepData{stepIdx+1} = avgTrace;
    end
end


% Ask where to save averaged result
saveDir = uigetdir(pwd, 'Select folder to save averaged .ma steps');

% Create a figure to hold all overlaid averaged traces
overlayFig = figure('Name', 'Averaged Traces', 'NumberTitle', 'off');
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
    
    % Copy a template raw .ma file (you must set this path to a raw file you know works)
templateFile = 'C:\Users\ah918\Desktop\VC E4_TR_DMSO\2025.06.03_000\slice_000\cell_001\OGS_ActInact_PNSub_000\000\Clamp1.ma';
destFile = fullfile(stepFolder, 'Clamp1.ma');

copyfile(templateFile, destFile);

% Overwrite /data dataset in copied file with averaged trace
fileId = H5F.open(destFile, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
dsetId = H5D.open(fileId, '/data');
H5D.write(dsetId, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', avgTrace);
H5D.close(dsetId);
H5F.close(fileId);

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


