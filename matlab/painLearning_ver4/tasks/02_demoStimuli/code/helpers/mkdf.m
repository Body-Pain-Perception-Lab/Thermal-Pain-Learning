function mkdf(vars,Results)

matName1 = strcat(vars.OutputFolder, vars.DataFileName(1:end-16), '_summary.mat');
csvName1 = strcat(vars.OutputFolder, vars.DataFileName(1:end-16), '_thresholds.csv');
csvName2 = strcat(vars.OutputFolder, vars.DataFileName(1:end-16), '_tcs.csv');

if ~exist(matName1) %#ok<EXIST>
    thresholds = [];
    tcsData = [];
    sots = [];
else
    load(matName1, 'thresholds', 'tcsData', 'sots')
end
    
trialCounter(:,1) = repmat(1:vars.nTrials,1,vars.nQstMeasures);
latest_thresholds = [Results.SubID  Results.trialN Results.thresholdN trialCounter Results.baseline Results.ThresholdT Results.ThresholdM];
latest_tcsData = Results.tcsData;
latest_sots = [Results.SOT_trial Results.SOT_jitter Results.SOT_stimOn Results.SOT_stimOff Results.SOT_ITI Results.TrialDuration];

thresholds = [thresholds; latest_thresholds];
tcsData = [tcsData; latest_tcsData];
sots = [sots; latest_sots];

% save summary mat file
save(matName1, 'thresholds', 'tcsData', 'sots');

% save as .csv
csvwrite(csvName1, thresholds);
csvwrite(csvName2, tcsData);
%struct2csv(Results, csvName);