function make_regression_contrasts(r_file_path, spm_file_path)
    % Load the R.mat file
    load(r_file_path, 'names', 'R');

    % Check if the required variables are present
    if ~exist('names', 'var') || ~exist('R', 'var')
        error('The loaded R.mat file does not contain the required "names" and "regressors" variables.');
    end

    regressors = R;
    % Load the SPM.mat file
    load(spm_file_path);

    % Ensure that the regressors in R.mat are in the SPM design matrix
    
    % note that in SPM multiple regression first col is always the
    % constant, so we need to shift everything by 1 
    
    if size(SPM.xX.X, 2)-1 < size(regressors, 2) || ~all(all(SPM.xX.X(:, 2:size(regressors, 2)+1) == regressors))
        error('The regressors in R.mat do not correspond to the beginning of the design matrix in SPM.mat.');
    end

    % Loop through the regressors and create both positive and negative contrasts
    njobs = 1;
    
    for i = 1:length(names)
        
        % positive contrast per regressor
        contrast_vector = zeros(1, size(SPM.xX.X, 2));
        
        contrast_vector(i+1) = 1;
        
        matlabbatch{njobs}.spm.stats.con.spmmat = {spm_file_path};
        matlabbatch{njobs}.spm.stats.con.consess{1}.tcon.name = sprintf('pos_%s',names{i});
        matlabbatch{njobs}.spm.stats.con.consess{1}.tcon.weights = contrast_vector;
        matlabbatch{njobs}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{njobs}.spm.stats.con.delete = 0;
    
        njobs = njobs+1;
        
        % put run job here
        % negative contrast per regressor
    
         % positive contrast per regressor
        contrast_vector = zeros(1, size(SPM.xX.X, 2));
        
        contrast_vector(i+1) = -1;
        
        matlabbatch{njobs}.spm.stats.con.spmmat = {spm_file_path};
        matlabbatch{njobs}.spm.stats.con.consess{1}.tcon.name = sprintf('neg_%s',names{i});
        matlabbatch{njobs}.spm.stats.con.consess{1}.tcon.weights = contrast_vector;
        matlabbatch{njobs}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{njobs}.spm.stats.con.delete = 0;
        
        njobs = njobs+1;
        
    end
    
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
