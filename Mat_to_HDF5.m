% Step 1: Select the top-level folder (e.g., VC E4_TR_DMSO)
rootDir = uigetdir(pwd, 'Select top-level folder (e.g., VC E4_TR_DMSO)');

if rootDir == 0
    disp('User cancelled.');
    return;
end

% Step 2: Find all Clamp1.ma files recursively
maFiles = dir(fullfile(rootDir, '**', 'Clamp1.ma'));

fprintf('Found %d Clamp1.ma files.\n', numel(maFiles));

% Step 3: Process each file
for k = 1:numel(maFiles)
    thisFile = fullfile(maFiles(k).folder, maFiles(k).name);
    fprintf('\nProcessing: %s\n', thisFile);

    % Step 3a: Try to identify if already HDF5
    try
        h5info(thisFile);
        fprintf('✅ Already HDF5. Skipped.\n');
        continue;  % Don't convert again
    catch
        % Not HDF5 — proceed to try loading as MAT file
    end

    % Step 3b: Try loading as a MATLAB .mat-style .ma file
    try
        s = load(thisFile);  % Will fail if not .mat format
        if isfield(s, 'data')
            data = s.data;
        else
            warning('❌ No "data" field found. Skipping.');
            continue;
        end
    catch ME
        warning('❌ Could not load file as .mat: %s\n', ME.message);
        continue;
    end

    % Step 3c: Delete original file and replace with HDF5
    try
        delete(thisFile);  % Remove old .mat version

        % Ensure data is a row vector (or modify if needed)
        data = data(:)';  % Force row vector

        h5create(thisFile, '/data', size(data));
        h5write(thisFile, '/data', data);
        fprintf('✅ Converted to HDF5.\n');
    catch ME
        warning('❌ Failed to convert %s: %s', thisFile, ME.message);
    end
end

disp('✨ Batch conversion complete.');
