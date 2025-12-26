function [EEG, T] = cleanTrial(EEG)
%% main cleaning function for pre-processing EEG data for DOD
% 1. remove ref channel
% 2. low pass filter
% 3. high pass filter
% 4. downsample
% 5. remove channels that is 3 std above the median rms of all channels
%%

% remove the flat ref channel so it is not listed as removed channel below
EEG = pop_select(EEG,'nochannel',{'5Z'});

% remove channels with flatlines more than 5seconds
EEG = clean_flatlines(EEG); 
removedChans_idx = find(strcmp({EEG.chaninfo.removedchans.type},'EEG'));
if ~isempty(removedChans_idx)
    T.chanNum_removedbyFlatline = length(removedChans_idx);
    T.channels_removedbyFlatline = {EEG.chaninfo.removedchans(removedChans_idx).labels};
else
    T.chanNum_removedbyFlatline = [];
    T.channels_removedbyFlatline = [];
end

% low pass
EEG = pop_eegfiltnew(EEG, [],85,112,0,[],0);   %lowpass cufoff at 100Hz with transition band of 30Hz and passband edge of 85Hz

% high pass
EEG = pop_eegfiltnew(EEG,'locutoff',2.5);    %highpass cutoff at 1.5Hz with transition band of 2 Hz and passband edge at 2.5Hz

% down sample from 1024 to 256
EEG = pop_resample(EEG, 256);

% % remove line noise
%EEG = cleanline(EEG,'newversion',1);
[EEG.data, EEG.etc.zapline.config, EEG.etc.zapline.analyticsResults] = clean_data_with_zapline_plus(EEG.data, EEG.srate,...
    'minfreq',45,'coarseFreqDetectPowerDiff',2.5,'plotResults',0);
[EEG, EEG.checkset_hist] = eeg_checkset(EEG);

% remove channels that is 3 std above the mean rms of the channel
rms_c = rms(EEG.data');
median_c = median(rms_c);
mean_c = mean(rms_c);
std_c = std(rms_c);
lower_bound = median_c - 3 * std_c;
upper_bound = median_c + 3 * std_c;
outRMS_chan = find(rms_c < lower_bound | rms_c > upper_bound);
rms_removed_ch = []; %tracker for RMS channel removal
channelRemoved1_num = [];
if ~isempty(outRMS_chan)
    rms_removed_ch = {EEG.chanlocs(outRMS_chan).labels};
    EEG = pop_select( EEG, 'rmchannel',rms_removed_ch);
    channelRemoved1_num = length(outRMS_chan);
else
    channelRemoved1_num = 0;
end

% % ICC %%%
% % apply iCanClean parameters
% params = [];
% params.rhoSqThres_source = .8;
% params.filtYtype = 'Notch';
% params.filtYfreq = [3 55];
% params.plotStatsOn = false;
% %params.rerefX = 'yes-fullrank';
% if ~contains(EEG.condition,'rest')
%     EEG = iCanClean(EEG, [1:length(EEG.chanlocs)], [1:length(EEG.chanlocs)] ,0, params);
% end
% 
% % remove channels using clean_channels
% rms_removed_chans = {EEG.chaninfo.removedchans(strcmp({EEG.chaninfo.removedchans.type},'EEG')).labels};
% [EEG,~] = clean_channels(EEG,.7,[],[],[],150,[]); % corr_threshold @ .7
% bad_chans = EEG.chaninfo.removedchans(strcmp({EEG.chaninfo.removedchans.type},'EEG'));
% if ~isempty(bad_chans)
%     cleanraw_removed_ch = {bad_chans.labels};
%     cleanraw_removed_ch = setdiff(cleanraw_removed_ch, rms_removed_chans, 'stable');
% else
%     cleanraw_removed_ch = [];
% end

% tracked removed channels
T.subject = EEG.subject;
T.condition = EEG.condition;
T.zaplineP = 0;
if isfield(EEG.etc,'zapline') && ~isempty(EEG.etc.zapline.analyticsResults.noisePeaks)
    T.zaplineUsed = 1;
    T.zaplineP = EEG.etc.zapline.analyticsResults.noisePeaks;
end
T.channels_kept = EEG.nbchan;
T.chanNum_removedbyRMS = channelRemoved1_num;
T.channels_removedbyRMS = rms_removed_ch;
%T.chanNum_removedbyCleanChan = length(cleanraw_removed_ch);
%T.channels_removedbyCleanChan = cleanraw_removed_ch;

end