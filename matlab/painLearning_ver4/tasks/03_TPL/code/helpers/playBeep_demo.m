
%%
device = PsychPortAudio('GetDevices',[],[]);
for i = 1:length(device)
    if strcmp(device(i).DeviceName, 'Speakers (Realtek(R) Audio)')
        deviceId = i;
    end
end
%%
clear; clc
InitializePsychSound(1); % initialize sound driver
% parameters
deviceIndex = 1;
nchannels = 2;
freq = 44100;
reps = 1;
beepLenghtSecs = 0.25;
startCue = 0; % start 0 = immediately 
waitForDeviceStart = 1;
pahandle = PsychPortAudio('Open',5);
%,deviceIndex,[],1,freq,nchannels);
%%
PsychPortAudio('Volume',pahandle,0.5); % set volume to 50%

myBeep = MakeBeep(1000,beepLengthSecs,freq); % generate beep
myBeep = repmat(myBeep,2,1); % to make it stereo

% Fill the audio playback buffer with the data, doubled for stereo
% presentation
PsychPortAudio('FillBuffer', pahandle, myBeep);

% Start Audio palyback
PsychPortAudio('Start', pahandle, reps, startCue, waitForDeviceStart)

% Wait for the beep to end
PsychPortAudio('Stop', pahandle, 1, 1);


PsychPortAudio('Close',pahandle);
