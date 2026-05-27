function sheet = makeSheetName(structName, fieldName)

    sheet = [structName '_' fieldName];
    sheet = regexprep(sheet, '[][\\/:*?]', '_');

    if strlength(sheet) > 31
        sheet = extractBefore(sheet, 32);
    end
end
