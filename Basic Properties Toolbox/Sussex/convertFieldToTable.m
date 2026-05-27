function T = convertFieldToTable(data, meta)

    n = numel(data);

    file_name   = repmat({meta.file_name},   n, 1);
    file_path   = repmat({meta.file_path},   n, 1);
    struct_name = repmat({meta.struct_name}, n, 1);
    field_name  = repmat({meta.field_name},  n, 1);
    index       = (1:n)';

    if isnumeric(data) || islogical(data)
        values = data(:);

    elseif iscell(data)
        values = data(:);

    elseif ischar(data) || isstring(data)
        values = repmat({data}, n, 1);

    else
        error('Unsupported data type: %s', class(data));
    end

    T = table( ...
        file_name, ...
        file_path, ...
        struct_name, ...
        field_name, ...
        index, ...
        values ...
    );
end
