function batchFixClampFilesInteractive()
    % Select folders and template
    parentFolder = uigetdir('', 'Select parent folder containing step folders');
    if parentFolder == 0, return; end

    referenceFile = uigetfile('*.ma', 'Select reference HDF5 file (raw trace with full metadata)');
    if referenceFile == 0, return; end
    refPath = fullfile(referenceFile);

    outputFolder = uigetdir('', 'Select folder to save fixed files');
    if outputFolder == 0, return; end

    % Step folders (e.g., 000, 001, ...)
    subdirs = dir(fullfile(parentFolder, '*'));
    subdirs = subdirs([subdirs.isdir] & ~startsWith({subdirs.name}, '.'));

    for k = 1:length(subdirs)
        try
            stepName = subdirs(k).name;
            inputPath = fullfile(parentFolder, stepName, 'Clamp1.ma');

            if ~isfile(inputPath)
                warning('File not found: %s', inputPath);
                continue;
            end

            % Confirm it’s HDF5 and readable
            try
               try
    fid = H5F.open(inputPath, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    dset = H5D.open(fid, '/data');
    avgTrace = H5D.read(dset);
    H5D.close(dset);
    H5F.close(fid);
catch
    warning('Skipping unreadable file: %s', inputPath);
    continue;
end

            catch
                warning('Skipping non-HDF5 file: %s', inputPath);
                continue;
            end

            % Create output file path
            outPath = fullfile(outputFolder, stepName, 'Clamp1.ma');
            if ~exist(fileparts(outPath), 'dir')
                mkdir(fileparts(outPath));
            end

            % Create new file and write /data
            if isfile(outPath)
                delete(outPath);  % avoid append issues
            end
            h5create(outPath, '/data', size(avgTrace), 'Datatype', class(avgTrace));
            h5write(outPath, '/data', avgTrace);

            % Copy metadata from reference
            copyHDF5Metadata(refPath, outPath);

        catch ME
            warning('Failed to process %s: %s', subdirs(k).name, ME.message);
        end
    end

    disp('Batch conversion complete.');
end

function copyHDF5Metadata(srcFile, destFile)
    % Read from /info and /MetaArray
    src = matfile(srcFile, 'Writable', false);
    try
        info = h5info(srcFile, '/info');
        copyGroup(srcFile, destFile, '/info');
    catch
        warning('No /info group found in %s', srcFile);
    end
    try
        meta = h5read(srcFile, '/MetaArray');
        h5create(destFile, '/MetaArray', size(meta), 'Datatype', class(meta));
        h5write(destFile, '/MetaArray', meta);
    catch
        warning('No /MetaArray dataset found in %s', srcFile);
    end
end

function copyGroup(srcFile, destFile, groupPath)
    info = h5info(srcFile, groupPath);

    for i = 1:length(info.Groups)
        srcGroup = info.Groups(i);
        newPath = srcGroup.Name;
        try
            h5createGroup(destFile, newPath);
        catch
            % group might already exist
        end
        copyGroup(srcFile, destFile, newPath);
    end

    for i = 1:length(info.Datasets)
        ds = info.Datasets(i);
        fullPath = fullfile(groupPath, ds.Name);
        data = h5read(srcFile, fullPath);
        h5create(destFile, fullPath, size(data), 'Datatype', class(data));
        h5write(destFile, fullPath, data);
    end
end

function h5createGroup(filename, groupname)
    % Create an empty group (must use low-level API)
    fid = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
    gid = H5G.create(fid, groupname, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    H5G.close(gid);
    H5F.close(fid);
end
