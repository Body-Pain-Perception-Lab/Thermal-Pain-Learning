function mpm_vbm(analysis_name,regressors, basepath)
% mpm_vbm(analysis_name,regressors)
% Function to set up stat analysis nad run mpm vbq
%
% input
%   analysis_name           name of the analysis, as string
%   regressors              full path to regressors file. mat file containing matrix 'R' and cell array 'names'
%

% FOR TESTING
% analysis_name = 'respiroceptAnalysis_1';
% regressors = '/mnt/fast_scratch/VBM/analysis/respiroception/regressors/respirocept_full_regressors_clean.mat';
if nargin < 3

    vbm_analysis_basepath =  '/mnt/fast_scratch/VBM/analysis/respiroception/';

else-24.00
    vbm_analysis_basepath = basepath;
end
% MPM data directory
hMRI_spatproc_folder = fullfile(filesep,'mnt','raid0','scratch','BIDS','derivatives','hMRI','spatproc','Results',filesep);

% Load regressors, mat file containing matrix 'R' and cell array 'names'
load(regressors)
subset = R(:,1);            % subset of sIDs who have data for this analysis



%% Load TIV file
tiv_file_path = fullfile(filesep,'mnt','raid0','scratch','BIDS','derivatives','summaries','mpm','TIV_vmp_all.csv'); % tab delimited to get sID
tivData = readtable(tiv_file_path);
tiv_ids = tivData(:,9);                 % extract sID in order
tivData = readtable(tiv_file_path,'Delimiter',{'\t',','},'MultipleDelimsAsOne',true);       % comma delimited to get TIV values
tivData = horzcat(tiv_ids,tivData);                                                         % join to sID in same order
tiv_ids = table2array(tiv_ids);                     
set = tiv_ids;
idxSubsetToSet  = findIdxSubject(set, subset);                                              % extract subset with task data on interest
tiv_ordered = tivData(idxSubsetToSet(1,:),:);                                       
tiv_ordered(:,2) = [];                                                                      % remove filepath column
tiv_ordered = table2array(tiv_ordered);   
tiv_ordered(:,end+1) = tiv_ordered(:,2) + tiv_ordered(:,3) + tiv_ordered(:,4);
tiv_mean = mean(tiv_ordered(:,end));
tiv_ordered(:,end) = tiv_ordered(:,end)-tiv_mean;

% Add TIV to R array
tiv_final = tiv_ordered(:,[1,end]);

%%

[matchingIDs, indices_R, indices_TIV] = intersect(R(:,1), tiv_final(:,1));
matched_R=R(indices_R, :);
matched_TIV = tiv_final(indices_TIV, 2:end); %Exclude the ID column from B as its redundant

R = [matched_R, matched_TIV];

%% R = horzcat(R, tiv_final(:,end)); this solution does not allow for cases where there is more behavior data than MRIs

names{end+1}='tiv';





%% Get desired wap1 images for all subjects
inputMapTypes =  {'R1', 'MT', 'R2s'};
time_string = datestr(now,'dd_mm_yy-HH_MM');



for mapLoop = 1:3
    
    inputMapType    = inputMapTypes{mapLoop};
    
    switch inputMapType
        case 'R1'
            
            % Create analysis directory & save workspace there
            analysis_name_R1 = fullfile([analysis_name, filesep, inputMapType,'_',time_string]);
            analysis_fullpath = fullfile(vbm_analysis_basepath, analysis_name_R1);
            analysis_dir_R1 =  {{analysis_fullpath}};
            
            filt= '^wap1.*\_R1.nii*';         
            % select map files
            [files_c1, ~] = spm_select('FPListRec', hMRI_spatproc_folder, filt);
            %files_c1 = [files_c1, repmat(',1',size(files_c1,1),1)];
            
            % Get sIDs from MPM maps, in order
            subjectNoList = [];
            for i = 1: size(files_c1,1)
                newStr      = extractAfter(files_c1(i,:),'SubjectNo');
                subjectNo   = newStr(1:4);
                subjectNoList = [subjectNoList; subjectNo];
            end
            subjectNoList_num = str2num(subjectNoList);
            
            % Select only map files for sIDs who have task/analysis data
            set = subjectNoList_num;
            idxSubsetToSet  = findIdxSubject(set, subset);
            
            maps_ordered    = files_c1(idxSubsetToSet(1,:),:);
            
            files_R1_cell = {cellstr(maps_ordered)};
            
        case 'MT'
                        
            % Create analysis directory
            analysis_name_MT = ([analysis_name, filesep, inputMapType,'_',time_string]);
            analysis_fullpath = fullfile(vbm_analysis_basepath, analysis_name_MT);
            analysis_dir_MT =  {{analysis_fullpath}};
            
            filt= '^wap1.*\_MTsat.nii*';
            [files_c1, ~] = spm_select('FPListRec', hMRI_spatproc_folder, filt);
            %files_c1 = [files_c1, repmat(',1',size(files_c1,1),1)];
            
            % Get sIDs from MPM maps, in order
            subjectNoList = [];
            for i = 1: size(files_c1,1)
                newStr      = extractAfter(files_c1(i,:),'SubjectNo');
                subjectNo   = newStr(1:4);
                subjectNoList = [subjectNoList; subjectNo];
            end
            subjectNoList_num = str2num(subjectNoList);
            
            % Select only map files for sIDs who have task/analysis data
            set = subjectNoList_num;
            idxSubsetToSet  = findIdxSubject(set, subset);
            maps_ordered    = files_c1(idxSubsetToSet(1,:),:);
            
            files_MT_cell = {cellstr(maps_ordered)};
            
        case 'R2s'
            
            % Create analysis directory
            analysis_name_R2s = ([analysis_name, filesep, inputMapType,'_',time_string]);
            analysis_fullpath = fullfile(vbm_analysis_basepath, analysis_name_R2s);
            analysis_dir_R2s =  {{analysis_fullpath}};
            
            
            filt= '^wap1.*\_R2s_OLS.nii*';
            [files_c1, ~] = spm_select('FPListRec', hMRI_spatproc_folder, filt);
            %files_c1 = [files_c1, repmat(',1',size(files_c1,1),1)];
            
            % Get sIDs from MPM maps, in order
            subjectNoList = [];
            for i = 1: size(files_c1,1)
                newStr      = extractAfter(files_c1(i,:),'SubjectNo');
                subjectNo   = newStr(1:4);
                subjectNoList = [subjectNoList; subjectNo];
            end
            subjectNoList_num = str2num(subjectNoList);
            
            % Select only map files for sIDs who have task/analysis data
            set = subjectNoList_num;
            idxSubsetToSet  = findIdxSubject(set, subset);
            maps_ordered    = files_c1(idxSubsetToSet(1,:),:);
            
            files_R2s_cell = {cellstr(maps_ordered)};
            
            finalBeh = subset(idxSubsetToSet(2,:));
    end
    
    if ~isfolder(analysis_fullpath)
        mkdir(analysis_fullpath)
    end
    save([analysis_fullpath,filesep,'workspace.mat'])
    
    clear files_c1 maps_ordered analysis_fullpath
end

% Clean up workspace
clear filt i idxSubsetToSet inputMapType inputMapTypes mapLoop newStr set subjectNo subjectNoList subjectNoList_num subset time_string


%% Save cleaned up names and R to regressors file

% exclude subjects from regressors without brain data

idxSubsetToSet = findIdxSubject(R(:,1), finalBeh);
R= R(idxSubsetToSet(1,:),:);
%% Remove sID from R and names
names(:,1) =  [];
R(:,1) =  [];

%% save final
[filepath,regfname] = fileparts(regressors);
regressors_final = fullfile(filepath,[regfname,'_tiv.mat']);
%if ~exist(regressors_final,'file')
    save(regressors_final,'R','names');
%end

regressors_final_cell2SPM =  {{regressors_final}};


%% Load job
create_design_mat_job
inputs{1}=matlabbatch;

%% Run
% initialise without graphics windows
spm('defaults', 'fmri')
%spm_jobman('initcfg')
spm_get_defaults('cmdline',true)
analysisError        = [];

try
    spm_jobman('run', inputs);
    inputs     = [];
catch ME
    analysisError = {analysisError, ME};
    
    rethrow(ME);
end

end




%% Helpers
% Helper Function
% function idxSubsetToSet = findIdxSubject(set, subset)
%     idxSubsetToSet = nan(2,length(subset));
%     count = 1;
%     for iter = 1:length(subset)
%         if ~isempty(find(set == subset(iter)))
%             idxSubsetToSet(1,count) = find(set == subset(iter));
%             idxSubsetToSet(2,count) = iter;
%             count = count + 1;
%         end
%     end
% 
%     % If any NaNs present, remove and reshape
%     if sum(isnan(idxSubsetToSet(:))) > 0
%         idxSubsetToSet(isnan(idxSubsetToSet)) = [];
%         howManyLeft = length(idxSubsetToSet)/2;
%         idxSubsetToSet = reshape(idxSubsetToSet, [2, howManyLeft]);
%     end
% end
% % 
% % function idxSubsetToSet = findIdxSubject(set, subset)
% %     idxSubsetToSet = []; % Start with an empty array
% %     for iter = 1:length(subset)
% %         idx = find(set == subset(iter), 1); % Find the index of the matching element in set
% %         if ~isempty(idx)
% %             idxSubsetToSet(end + 1, 1) = idx; % Add the matching index from set
% %             idxSubsetToSet(end, 2) = iter;    % Add the current index from subset
% %         end
% %     end
% % end
% 
% function idxSubsetToSet = findIdxSubjectStr(set, subset)
%     idxSubsetToSet = []; % Start with an empty array
%     for iter = 1:length(subset)
%         subsetStr = sprintf('%04d', subset(iter)); % Convert subset value to string with leading zeros
%         idx = find(strcmp(set, subsetStr), 1); % Find the index of the matching element in set
%         if ~isempty(idx)
%             idxSubsetToSet(end + 1, 1) = idx; % Add the matching index from set
%             idxSubsetToSet(end, 2) = iter;    % Add the current index from subset
%         end
%     end
% end

% function idxSubsetToSet = findIdxSubject(set, subset)
%     idxSubsetToSet = nan(2,length(subset));
%     count = 1;
%     for iter = 1:length(subset)
%         idx = find(set == subset(iter), 1);
%         if ~isempty(idx)
%             idxSubsetToSet(1,count) = idx;
%             idxSubsetToSet(2,count) = iter;
%             count = count + 1;
%         end
%     end
%     
%     % Remove NaN columns
%     idxSubsetToSet(:, isnan(idxSubsetToSet(1,:))) = [];
% end
function idxSubsetToSet = findIdxSubject(set, subset)
    idxSubsetToSet = nan(2,length(subset));
    count = 1;
    for iter = 1:length(subset)
        idx = find(set == subset(iter), 1);
        if ~isempty(idx)
            idxSubsetToSet(1,count) = idx;
            idxSubsetToSet(2,count) = iter;
            count = count + 1;
        end
    end
    
    % Remove NaN columns
    idxSubsetToSet(:, isnan(idxSubsetToSet(1,:))) = [];
end
