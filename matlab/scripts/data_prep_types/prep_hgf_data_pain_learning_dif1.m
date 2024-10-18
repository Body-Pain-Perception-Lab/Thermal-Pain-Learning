function prep_hgf_data_pain_learning_dif1(nSubjects, removers)
%prep_hgf_data(nSubjects)
%   Prepare pain learning task data for several subjects for the HGF workshop



%% setup
if nargin < 1
    nSubjects = 10;
end

dataFolder      = fullfile('C:\Users\Jespe\OneDrive\Skrivebord\pain-learning\TPL\TPL');              % <--- EDIT THE DATA FOLDER HERE
% subDirectories  = dir(dataFolder);
FileList  = dir(dataFolder);
% searchString    = 'TPL_v*-1_0';                                                                     % <--- AND THE FILE NAME
dataFolderOut   = fullfile('C:\Users\Jespe\OneDrive\Skrivebord\pain-learning\TPL\realdeal');                              % <--- EDIT OUTPUT FOLDER

tplSequenceFolder = fullfile('C:\Users\Jespe\OneDrive\Skrivebord\pain-learning\bachelor\pain experiment script\painLearning_ver4\tasks\03_TPL');
addpath (tplSequenceFolder)

a = 0
u=[];
u_cues = [];
u_densProb = [];
y=[];
y_rts=[];
y_conf=[];
y_conf_rts=[];
stim1 = [];
stormdb1 = [];
y_pred = [];
predacc1 = [];
v = [];
%% Loop over subject

for subLoop = 3 : nSubjects+2 %skip dummy folders     

        resultsFile = fullfile(FileList(subLoop).folder, FileList(subLoop).name);
        load(resultsFile)
        
        %% Get subject name and load cuo0Prediction and cueProbabilityOutput - These script excerpts should reflect TPL_loadParams.m 
        
        split = strsplit(FileList(subLoop).name,{'_','r.mat'},'CollapseDelimiters',true);
        vars.subNo = str2double(split{3});
        
        test = str2double(split(:,3)) == removers
        
        if (sum(test)==1)
            a = a+1
           continue 
        end
        v = [v, str2double(split{2}(2))];
%         if(v(subLoop-2) == 1)
%             continue
%         end
        
        if mod(vars.subNo,2)
            chosenSequence = 'sequence1_tgi.mat';
            vars.seqN = 1;
        else
            chosenSequence = 'sequence4_tgi.mat';
            vars.seqN = 4;
        end
        [Output] = load(chosenSequence);
        
        vars.conditionSequence   = Output.cueProbabilityOutput(:,2);     % 1 cue_0 valid
                                                                         % 2 cue_1 valid
                                                                         % 3 cue_0 invalid    
                                                                         % 4 cue_1 invalid      
                                                                         % 5 non-predictive
                                                                         % 6 cue_0 TGI
                                                                         % 7 cue_1 TGI

        vars.cue0Prediction      = Output.cueProbabilityOutput(:,12);           %0 NP, 1 cue_0->Warm, 2 cue_0->Cold% 
        vars.cueSequence         = Output.cueProbabilityOutput(:,5);            % sequence of cues [0|1]
        vars.TempQual            = Output.cueProbabilityOutput(:,10);           % temp outcome. 0 Cold  or 1 Warm
        vars.trialSequence       = Output.cueProbabilityOutput(:,6);            %trial type          [1, 2, 3]           % Valid / Invalid / TGI
        vars.desiredProb         = Output.cueProbabilityOutput(:,7);
        
        %% check results for missing data
%         missing_data = unique([find(1-Results.vasYN); find(isnan(Results.predResp))]);
%         missing_data = or(isnan(Results.SOT_trial), isnan(Results.predResp));
%         sum(missing_data)
%         %% clean up
%         
%         Results.predRT(missing_data)         = [];
%         
%         vars.conditionSequence(missing_data)= [];
%         vars.cue0Prediction(missing_data)   = [];
%         vars.cueSequence(missing_data)      = [];
%         vars.TempQual(missing_data)         = [];
%         vars.trialSequence(missing_data)    = [];
%         vars.desiredProb(missing_data)      = [];
%         
%         Results.predAcc(missing_data)       = [];
%         Results.predResp(missing_data)      = [];
%         Results.vasRT(missing_data,:)       = [];
%         Results.vasResp(missing_data,:)     = [];

        %% get variables
        
        missed_trials       = Results.predRT >= 3; %% missed trials vector
        nonpred_trials      = vars.cue0Prediction == 0; % 1 non-predicitve trial
        outcomes            = vars.TempQual;   % 0 cold, 1 warm, 2 TGI
        cues                = vars.cueSequence;  % 0 cold, 1 warm, 2 TGI
        responses           = Results.predResp;  %% 1 = warm, 0 = cold
        cue_valid           = vars.trialSequence; % 1 valid, 2 invalid, 3 tgi
        predacc             = Results.predAcc;
        
        response_rts        = Results.predRT;
        confidence_rts      = Results.vasRT;
        confidence_ratings  = Results.vasResp;
        stim                = Results.stim;
        stormdb             = split(:,3);
        %Check
        conditions          = vars.conditionSequence;   
        
        tgi_trials = or(conditions == 6, conditions == 7);
        des_prob = vars.desiredProb;
        
        for n = 1:length(des_prob)
            
            if tgi_trials(n)
                
                des_prob(n) = des_prob(n-1);
                
            end
            
        end
        %% get outcomes so that 1 = cueA -> cold and cueB -> warm, 2 = cueA -> warm and cueB -> cold
outcomes1 = outcomes;
dif = nan(length(outcomes1), 1);
        for i = 1:306
        if outcomes1(i) == 2;
            dif(i) = confidence_ratings(i,1)-confidence_ratings(i,2);
        end
        end
        
        for i = 1:306;
        if dif(i) > 0;
            outcomes1(i) = 0;
        elseif dif(i) < 0;
            outcomes1(i) = 1;
        end
        end
        
        cond1 = or(cues == 0 & outcomes1 == 0, cues == 1 & outcomes1 == 1);
        %cond1(cond1(outcomes ~= 2)) = NaN;
        
        cond2 = or(cues == 0 & outcomes1 == 1, cues == 1 & outcomes1 == 0);
        %cond2(cond2(outcomes ~= 2)) = NaN;
        
        inputs = nan(length(cond1), 1);
        inputs(cond1) = 0;
        inputs(cond2) = 1;
        inputs(outcomes1 == 2) = 0.5; %TGI
        
        %fills array with nan i different number of trials
        if length (inputs)<306
            inputs(length(inputs):306)=nan;
            des_prob(length(des_prob):306)=nan;
        end
        
        u = [u, inputs];
        u_densProb = [u_densProb, des_prob];
        
        %% get responses so that 1 = cueA -> cold and cueB -> warm, 2 = cueA -> warm and cueB -> cold, 3=TGI 
        cond1 = or(cues == 0 & responses == 0, cues == 1 & responses == 1);
        %cond1(cond1(outcomes ~= 2)) = NaN;
        
        cond2 = or(cues == 0 & responses == 1, cues == 1 & responses == 0);
        %cond2(cond2(outcomes ~= 2)) = NaN;
        
        cond3 = responses==2; %TGI
        
        inputs = nan(length(cond1), 1);
        inputs(cond1) = 0;
        inputs(cond2) = 1;
%         inputs(outcomes == 2) = NaN; % NoTGI
%         inputs_shift = circshift (inputs,1); %TGI
%         inputs (cond3) = inputs_shift(cond3); %TGI
%         
        %fills array with nan i different number of trials
%         if length (inputs)<306
%             inputs(length(inputs):306)=nan;
%             response_rts(length(response_rts):306)=nan;
%             confidence_ratings(length(confidence_ratings):306,:)=nan;
%             confidence_rts(length(confidence_rts):306,:)=nan;
%             responses(length(responses):306)=nan;
%             predacc(length(responses):306)=nan;
%             cues(length(cues):306)=nan;
%         end
        
        y = [y, inputs];
        predacc1 = [predacc1, predacc];
        y_pred = [y_pred, responses];
        y_rts = [y_rts, response_rts];
        y_conf (:,:,subLoop-2-a) = confidence_ratings;
        y_conf_rts (:,:,subLoop-2-a)= confidence_rts;
        stim1 = [stim1, stim];
        stormdb1 = [stormdb1, str2double(stormdb)];
        u_cues = [u_cues, cues];  %This was added due to the habituation effect, included on the perception model
        clear Results vars

        
%     end
end
v = (v-1)'
stormdb1 = stormdb1'
rmpath(tplSequenceFolder)

%% save as hgf-ready datafilea
thisSub         = 'multsub_TGI';
DataFileName    = ['HGFwksh_TPL_', thisSub,'.mat'];
save(strcat(dataFolderOut, filesep, DataFileName),  'u', 'u_cues','u_densProb','y', 'y_rts', 'y_conf', 'y_conf_rts', 'stim1', 'stormdb1','y_pred','predacc1','v');
disp(['Data saved as: ', DataFileName]);

