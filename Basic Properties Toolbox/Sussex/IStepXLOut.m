%% Select MAT-files via dialog (multi-select)
[fileNames, filePath] = uigetfile( ...
    '*.mat', ...
    'Select one or more MAT-files', ...
    'MultiSelect', 'on' ...
);

if isequal(fileNames, 0)
    error('No files selected.');
end

% Ensure consistent cell array
if ischar(fileNames)
    fileNames = {fileNames};
end

%% Output Excel file
outputXLS = fullfile(filePath, 'exported_data.xlsx');

%% Build file list (dir()-like structure)
files = struct([]);
for i = 1:numel(fileNames)
    files(i).name   = fileNames{i};
    files(i).folder = filePath;
end

%% Loop through each selected MAT-file
for f = 1:numel(files)

    matPath = fullfile(files(f).folder, files(f).name);
    data    = load(matPath);

    structNames = fieldnames(data);

    %% Loop through structures
    for s = 1:numel(structNames)

        structName = structNames{s};
        S = data.(structName);

        if ~isstruct(S)
            warning('%s in %s is not a struct — skipping', ...
                structName, files(f).name);
            continue
        end

        fieldNames = fieldnames(S);

        %% Loop through fields
        for k = 1:numel(fieldNames)

            fieldName = fieldNames{k};
            fieldData = S.(fieldName);

            % Traceability metadata
            meta.file_name   = files(f).name;
            meta.file_path   = matPath;
            meta.struct_name = structName;
            meta.field_name  = fieldName;

            try
                T = convertFieldToTable(fieldData, meta);
            catch ME
                warning('Failed on %s.%s (%s): %s', ...
                    structName, fieldName, files(f).name, ME.message);
                continue
            end

            sheetName = makeSheetName(structName, fieldName);
            writetable(T, outputXLS, 'Sheet', sheetName);

        end
    end
end

disp('Export complete.');
