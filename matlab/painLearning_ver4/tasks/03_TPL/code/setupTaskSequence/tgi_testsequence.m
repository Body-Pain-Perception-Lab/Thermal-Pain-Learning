% Create TGI stimuli specifications 
%% Load sequence to modify
sequenceN = 4;
load(['sequence', num2str(sequenceN), '.mat'])

%% Find how many blocks
blocksI = [find(cueProbabilityOutput(:,15) == 1); length(cueProbabilityOutput)];
blocksN = length(blocksI);
blocks = [];
for i = 2:blocksN
    block_temp = ones(blocksI(i)-blocksI(i-1),1)*i-1;
    blocks = [blocks;block_temp]; %#ok<AGROW>
end
blocks(length(cueProbabilityOutput)) = blocksN-1;
%% Update block numbers
cueProbabilityOutput(:,16) = blocks;

%% Block sizes
blocksSize = zeros(1,blocksN-1);
for i = 1:blocksN-1
    blocksSize(i) = sum(blocks==i); % tot number of stimuli in a block
end

%%
%stim_n           = blocksSize; 
block_type       = cueProbabilityOutput(blocksI(1:end-1),3); % block type
face_gender      = zeros(blocksN-1,1);
cue              = zeros(blocksN-1,1);
trial_type       = ones(blocksN-1,1)*3;
desired_prob     = cueProbabilityOutput(blocksI(1:end-1),7);  
effective_prob   = cueProbabilityOutput(blocksI(1:end-1),8);
block_volatility = cueProbabilityOutput(blocksI(1:end-1),9);
tgi_n_20         = round(blocksSize/100*15); % 20% of trials in a block
%tgi_n_even       = tgi_n_20 - mod(tgi_n_20,2); % even number of tgi trials per block
%tgi_n_even       = [8 2 4 4 4 4 8 4 4 2 4]; 
tgi_n_even       = [6 2 4 4 4 2 6 4 4 2 4]; 
tgi_1            = round(blocksSize/100*10); % after 10% of trials in a block 
blocksSizeTgi    = blocksSize + tgi_n_even;
%% define trial index (from 1 to total number of stimuli)
tgi_trial = zeros(blocksN-1,1);
tgi_trial_start = zeros(1,blocksN-1);
tgi_trial_end = zeros(1,blocksN-1);

for i = 1:blocksN-1
    if i==1
        tgi_trial_start(i) = blocksI(i)+tgi_1(i);
        tgi_trial_end(i) = blocksSize(i) + tgi_n_even(i);
    else
        tgi_trial_start(i) = blocksI(i)+tgi_1(i)+tgi_n_even(i-1);
        tgi_trial_end(i) = sum(blocksSize(1:i)) + sum( tgi_n_even(1:i));
    end
end
if sum(blocksSizeTgi) ~= tgi_trial_end(end)
    warning("tgi trial start/end may be defined incorrectly")
    return
end
%%
tgiI_trial = [];
for i = 1:blocksN-1
    tgiI_singleBlock    = randsample(tgi_trial_start(i):tgi_trial_end(i), tgi_n_even(i));
    tgiI_singleBlock    = sort(tgiI_singleBlock)';
    tgiI_trial       = [tgiI_trial; tgiI_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,1) = tgiI_trial;
%% define condition (6 or 7)
tgi_condition = [];
for i = 1:blocksN-1
    tgiC_singleBlock = [ones(tgi_n_even(i)/2,1)*6; ones(tgi_n_even(i)/2,1)*7];
    tgiC_singleBlock = Shuffle(tgiC_singleBlock);
    tgi_condition = [tgi_condition; tgiC_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,2) = tgi_condition;
%% define block type [1, 2, 3] 
tgi_blockType = [];
blockType = cueProbabilityOutput(blocksI(1:end-1),3);
for i = 1:blocksN-1
    tgiT_singleBlock = repelem(blockType(i),tgi_n_even(i),1);
    tgi_blockType    = [tgi_blockType; tgiT_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,3) = tgi_blockType;
%% face gender [0, 1, 2 = Not a face]
tgi_faceGender = ones(length(tgiI_trial),1)*2;
tgi_sequence(:,4) = tgi_faceGender;
%% cue [trial type 6 = cue 0, trial type 7 = cue 1]
tgi_cue = tgi_condition==7;
tgi_sequence(:,5) = tgi_cue;
%% trial validity [1 = valid, 2 = invalid, 3 = undefined]
tgi_validity =  ones(length(tgiI_trial),1)*3;
tgi_sequence(:,6) = tgi_validity;
%% desired prob
tgi_desProb = [];
desProb = cueProbabilityOutput(blocksI(1:end-1),7);
for i = 1:blocksN-1
    tgiD_singleBlock = repelem(desProb(i),tgi_n_even(i),1);
    tgi_desProb    = [tgi_desProb; tgiD_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,7) = tgi_desProb;

%% effective prob
tgi_effProb = [];
effProb = cueProbabilityOutput(blocksI(1:end-1),8);
for i = 1:blocksN-1
    tgiE_singleBlock = repelem(desProb(i),tgi_n_even(i),1);
    tgi_effProb    = [tgi_effProb; tgiE_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,8) = tgi_effProb;
%% block volatility [0 = stable, 1 = volatile]
tgi_blockVolatility = [];
blockVol = cueProbabilityOutput(blocksI(1:end-1),9);
for i = 1:blocksN-1
    tgiV_singleBlock = repelem(blockVol(i),tgi_n_even(i),1);
    tgi_blockVolatility = [tgi_blockVolatility; tgiV_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,9) = tgi_blockVolatility;
%% outcome [0 = cold, 1 = warm, 2 = both]
tgi_outcome = [];
outcome = 2;
for i = 1:blocksN-1
    tgiO_singleBlock = repelem(outcome,tgi_n_even(i),1);
    tgi_outcome    = [tgi_outcome; tgiO_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,10) = tgi_outcome;
%% predictive/non-predictive trial
tgi_predictive = [];
predictive = 2;
for i = 1:blocksN-1
    tgiP_singleBlock = repelem(predictive,tgi_n_even(i),1);
    tgi_predictive    = [tgi_predictive; tgiP_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,11) = tgi_predictive;

%% cue0PredictionSequence [0 = NP, 1 = cue_0->Warm, 2 = cue_0->Cold]
tgi_cue0Prediction = [];
blockCue0 = cueProbabilityOutput(blocksI(1:end-1),12);
for i = 1:blocksN-1
    tgiC0_singleBlock = repelem(blockCue0(i),tgi_n_even(i),1);
    tgi_cue0Prediction = [tgi_cue0Prediction; tgiC0_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,12) = tgi_cue0Prediction;

%% predictionTrialNext [0 = NP, 1 = P, 2 = Undefined] % 1 if there is a prediction trial after the ITI of this trial
tgi_predictiveNext = [];
predictiveNext = 2;
for i = 1:blocksN-1
    tgiPN_singleBlock = repelem(predictiveNext,tgi_n_even(i),1);
    tgi_predictiveNext = [tgi_predictiveNext; tgiPN_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,13) = tgi_predictiveNext;

%% reversalBlocksSequence; 
tgi_reversalBlocksSequence = [];
reversalBlocksSequence = 2;
for i = 1:blocksN-1
    tgiRBS_singleBlock = repelem(reversalBlocksSequence,tgi_n_even(i),1);
    tgi_reversalBlocksSequence = [tgi_reversalBlocksSequence; tgiRBS_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,14) = tgi_reversalBlocksSequence;
%% blockwiseTrialNumber                    % trial # within a block
tgi_sequence(:,15) = 1001:1000+length(tgiI_trial);
%% blockNumber                             % block #
tgi_blockNumber = [];
blockNumber = cueProbabilityOutput(blocksI(1:end-1),16);
for i = 1:blocksN-1
    tgiBN_singleBlock = repelem(blockNumber(i),tgi_n_even(i),1);
    tgi_blockNumber = [tgi_blockNumber; tgiBN_singleBlock];  %#ok<AGROW>
end
tgi_sequence(:,16) = tgi_blockNumber;

%% Redefine trial number of non-TGI stimuli
tgi_cueProbabilityOutput = cueProbabilityOutput;% zeros(size(tgi_sequence,1),size(tgi_sequence,2))];
%tgi_cueProbabilityOutput(:,1) = 1:length(tgi_cueProbabilityOutput);
% for s = 1:length(tgi_cueProbabilityOutput)
%     if any(tgi_sequence(:,1) == tgi_cueProbabilityOutput(s,1))
%         %index = find(tgi_sequence(:,1) == tgi_cueProbabilityOutput(s,1));
%         tgi_cueProbabilityOutput(s:end,1) = tgi_cueProbabilityOutput(s:end,1) + 1;
%     end
% end

%% Add tgi stimuli
cueProbabilityOutput = sortrows([tgi_cueProbabilityOutput; tgi_sequence]);
cueProbabilityOutput(:,1) = 1:length(cueProbabilityOutput);
if isempty(find([1;diff(cueProbabilityOutput(:,1))]~=1, 1))
    disp("Sequence has consecutive number trials")
else
    warning("Sequence does NOT have consecutive number trials")
    return
end

%% Check number of TGI stimuli per block
a = [blockCue0 desProb tgi_n_even'];
if sum(a(a(:,1)==0,3)) ~= sum(a(a(:,1)==1,3)) 
    warning("The number of TGI stimuli is different across blocks 0 and 1")
    return
elseif sum(a(a(:,1)==0,3)) ~= sum(a(a(:,1)==2,3)) 
    warning("The number of TGI stimuli is different across blocks 0 and 2")
    return
end

disp(['N TGI stimuli in block 0: ' num2str(sum(a(a(:,1)==0,3)))]);
disp(['N TGI stimuli in block 1: ' num2str(sum(a(a(:,1)==1,3)))]);
disp(['N TGI stimuli in block 2: ' num2str(sum(a(a(:,1)==2,3)))]);

%% Add whether a trial should be followed by ratings
cond_1 = find(cueProbabilityOutput(:,2) == 1);
cond_2 = find(cueProbabilityOutput(:,2) == 2);
cond_3 = find(cueProbabilityOutput(:,2) == 3);
cond_4 = find(cueProbabilityOutput(:,2) == 4);
cond_5 = find(cueProbabilityOutput(:,2) == 5);
cond_6 = find(cueProbabilityOutput(:,2) == 6);
cond_7 = find(cueProbabilityOutput(:,2) == 7);
condSize = length(cond_6);

cond_1r = randsample(cond_1,condSize);
cond_2r = randsample(cond_2,condSize);
cond_3r = randsample(cond_3,length(cond_3));
cond_4r = randsample(cond_4,length(cond_4));
cond_5r = randsample(cond_5,condSize);

conds = sort([cond_1r; cond_2r; cond_3r; cond_4r; cond_5r; cond_6; cond_7]);
ratings_yn = zeros(length(cueProbabilityOutput(:,1)),1);
for c = 1:length(ratings_yn)
    if any(cueProbabilityOutput(c,1)== conds)
        ratings_yn(c) = 1;
    end
end
cueProbabilityOutput(:,17) = ratings_yn; 
%ratings_yn(1:10)
%% Plot trials
%figure; hold on
%plot(cueProbabilityOutput(:,1),cueProbabilityOutput(:,10),'o')
%plot(cueProbabilityOutput(:,1),cueProbabilityOutput(:,17),'*')
%%
breaks = [100;220];
%save(['sequence', num2str(sequenceN) ,'_tgi.mat'],'blockParams','breaks','cueProbabilityOutput')