%model comparison
function c = get_matlab_data(mainScriptPath)

outcome_dir = fullfile(mainScriptPath,"matlab data");

BIDS_dir = fullfile(mainScriptPath,"csv_files2","Data");
subjectDirs = dir(BIDS_dir);

subjectDirs = subjectDirs([subjectDirs.isdir]); % Keep only directories

% Filter out subdirectories that start with "sub-"
subjectDirs = subjectDirs(startsWith({subjectDirs.name}, 'sub-'));


% Loop over each valid subject directory
for i = 1:length(subjectDirs)
      
    subjectName = subjectDirs(i).name;
    subjectDir = fullfile(BIDS_dir, subjectName);
    
    
    filePattern = '*tpl_beh.tsv';
    matchingFiles = dir(fullfile(subjectDir, filePattern));
    
    try
        tsvFilePath = fullfile(subjectDir, matchingFiles(1).name);
    catch
        continue
    end
    tableData = readtable(tsvFilePath, 'Delimiter', '\t', 'FileType', 'text');
    Results = table2struct(tableData, 'ToScalar', true);
     if ~iscell(Results.vasRT_2)
 
         string_column = arrayfun(@num2str, Results.vasRT_2, 'UniformOutput', false);
         string_column(isnan(Results.vasRT_2)) = {'NA'};
         Results.vasRT_2 = string_column
 
         string_column = arrayfun(@num2str, Results.vasRT_3, 'UniformOutput', false);
         string_column(isnan(Results.vasRT_3)) = {'NA'};
         Results.vasRT_3 = string_column
         
         
         string_column = arrayfun(@num2str, Results.vasResp_2, 'UniformOutput', false);
         string_column(isnan(Results.vasResp_2)) = {'NA'};
         Results.vasResp_2 = string_column
         
         string_column = arrayfun(@num2str, Results.vasResp_3, 'UniformOutput', false);
         string_column(isnan(Results.vasResp_3)) = {'NA'};
         Results.vasResp_3 = string_column
     end
    
    Results.vasRT=[Results.vasRT_1,Results.vasRT_2,Results.vasRT_3];
    Results.vasResp=[Results.vasResp_1,Results.vasResp_2,Results.vasResp_3];
    Results.targetT =[Results.targetT_1,Results.targetT_2,Results.targetT_3,Results.targetT_4,Results.targetT_5];


    digits = regexp(subjectName, '\d+', 'match');
    %needs the same name.
    outcome_sub = "TPL_v2-3_"+digits{1}+"r.mat"
    save(fullfile(outcome_dir,outcome_sub),"Results"); 

end


end