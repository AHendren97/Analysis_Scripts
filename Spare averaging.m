

% Script: average_pn_traces_auto_steps.m
% Description: Select first step file for each recording. Automatically loads
% subsequent steps from the same folder. Averages and saves each in
% 'OGS_Actinact_Pnsub_Avr/000', '001', etc.

% ----------- Parameters -----------
outputRoot = 'OGS_Actinact_Pnsub_Avr';  % Where to save output
Nsteps = 41;       % input('How many steps per P/N recording? ');
Nrecordings = input('How many recordings do you want to process? ');
ext = '.ma';  % Extension used by your files
% ----------------------------------

% Create output root if needed
if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

% Track last used path
lastPath = pwd;

for recIdx = 1:Nrecordings
    fprintf('\n--- Recording %d of %d ---\n', recIdx, Nrecordings);

    % Select only the first file
    [firstFile, subpath] = uigetfile(['*', ext], ...
        sprintf('Select FIRST step file for Recording %d', recIdx), lastPath);

    % Cancel check
    if isequal(firstFile, 0)
        fprintf('Selection cancelled for recording %d. Skipping.\n', recIdx);
        continue;
    end

    lastPath = subpath;  % Update for next recording
    [~, baseName, ~] = fileparts(firstFile);

    % Extract common prefix from filename (e.g. "cellA_step_000" → "cellA_step_")
    % This assumes format like: name_step_000.ma
    basePrefix = regexprep(baseName, '\d{3}$', '');  % strip trailing digits

    % Loop through each step file
    for stepIdx = 0:Nsteps-1
        stepStr = sprintf('%03d', stepIdx);
        fullName = fullfile(subpath, [basePrefix, stepStr, ext]);

        if ~isfile(fullName)
            fprintf('⚠️ Step file not found: %s (skipping)\n', fullName);
            continue;
        end

        % Load data
        sub = ephysIO(fullName);
        Time = sub.array(:, 1);
        Trace = sub.array(:, 2);

        % Collect for averaging
        if stepIdx == 0
            allTraces = zeros(length(Trace), Nsteps);
        end

        allTraces(:, stepIdx + 1) = Trace;
    end

    % Average across steps
    pn_avg = mean(allTraces(:, 1:Nsteps), 2, 'omitnan');

    % Save averaged file into the same folder as the input step files
saveName = fullfile(subpath, sprintf('pn_avg_R%02d%s', recIdx, ext));

% Save .ma file (MAT format renamed)
save([saveName(1:end-3), '.mat'], 'pn_avg', 'Time');
movefile([saveName(1:end-3), '.mat'], saveName);

fprintf('✓ Saved averaged trace to: %s\n', saveName);


    fprintf('✓ Saved averaged steps for Recording %d\n', recIdx);
end

fprintf('\n✅ All done! Files saved in %s\n', outputRoot);
