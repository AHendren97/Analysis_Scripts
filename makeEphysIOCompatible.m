function makeEphysIOCompatible(avgTrace, referenceFile, outputFile)
% Create a valid EphysIO-compatible .ma HDF5 file from averaged trace data.
%
% avgTrace     - the Nx1 or 1xN trace matrix
% referenceFile - path to a working .ma file to clone metadata from
% outputFile    - where to write the new .ma file

    if ~endsWith(outputFile, '.ma')
        outputFile = [outputFile '.ma'];
    end

    % Create base file with trace
    h5create(outputFile, '/data', size(avgTrace), 'Datatype', 'double');
    h5write(outputFile, '/data', avgTrace);

    % Copy metadata from the reference .ma file
    copyInfoFromReference(referenceFile, outputFile);

    disp(['Created compatible file: ' outputFile]);
end

function copyInfoFromReference(refFile, targetFile)
    info = h5info(refFile, '/info');
    recursivelyCopyGroup('/info', info, refFile, targetFile);
end

function recursivelyCopyGroup(groupPath, groupInfo, refFile, targetFile)
    % Recreate group in target file
    try
        h5createGroup(targetFile, groupPath);
    catch
        % Group might already exist (especially root '/info')
    end

    % Copy attributes
    for a = groupInfo.Attributes
        h5writeatt(targetFile, groupPath, a.Name, h5readatt(refFile, groupPath, a.Name));
    end

    % Copy datasets
    for d = groupInfo.Datasets
        datasetPath = [groupPath '/' d.Name];
        data = h5read(refFile, datasetPath);
        h5create(targetFile, datasetPath, size(data), 'Datatype', d.Datatype.Type);
        h5write(targetFile, datasetPath, data);

        % Copy dataset attributes
        for a = d.Attributes
            h5writeatt(targetFile, datasetPath, a.Name, h5readatt(refFile, datasetPath, a.Name));
        end
    end

    % Recursively copy subgroups
    for g = groupInfo.Groups
        recursivelyCopyGroup([groupPath '/' g.Name], g, refFile, targetFile);
    end
end

function h5createGroup(filename, groupPath)
    fid = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
    gid = H5G.create(fid, groupPath, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
    H5G.close(gid);
    H5F.close(fid);
end
