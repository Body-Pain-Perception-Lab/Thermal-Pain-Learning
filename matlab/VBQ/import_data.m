function dataTable = import_data(filePath)
    % Check if the file exists
    if exist(filePath, 'file') == 2
        % Read the CSV file and convert it to a table
        dataTable = readtable(filePath);
    else
        error('The specified file does not exist.');
    end
end
