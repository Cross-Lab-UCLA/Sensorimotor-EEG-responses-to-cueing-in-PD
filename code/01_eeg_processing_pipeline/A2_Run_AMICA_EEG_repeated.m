%% A2_Run_AMICA_EEG_repeated.m
%
%   This script runs AMICA ICA decomposition on merged EEG datasets for cluster-level
%   analysis. Bad segments are removed using joint probability and
%   amplitude criteria. Then, data is interpolated and re-referenced before
%   ICA. AMICA is performed on cleaned proxy data, and ICA weights are
%   added to the original data at the end.
%
%   1. Prompts user to select merged EEG datasets folder.
%   2. Loads each subject's EEG set and extracts walking/standing segments.
%   3. Removes bad segments using joint probability and amplitude thresholds, adjusting thresholds to retain 90% data.
%   4. Interpolates bad channels, adds reference channel, and re-references to average.
%   5. Runs AMICA ICA decomposition on cleaned data.
%   6. Adds ICA weights to original data and saves both cleaned and ICA-processed datasets.
%
% Outputs:
%   - AMICA-processed EEG datasets with ICA weights
%
% LM 122325
%% set up folders
clear all; clc;
eeglab; close all;
if ispc
    mainDir = 'E:\clab\DoD-Gait';
elseif ismac
    mainDir = '/Users/Leo/Git/DoD-Gait';
else
    mainDir = '/home/leo/Documents/DoD-Gait';
end
addpath(fullfile(mainDir,'code','func'))
disp('Load in 01_combined folder.')
workDir = uigetdir(fullfile(mainDir,'data'),'Load in 01_combined folder.');
label       = extractAfter(workDir,'combined_');
eegFiles = dir(fullfile(workDir,'*.set'));

%% main loop
for re = 1:20
    saveFolder =  fullfile(mainDir,'data',['02_AMICAed_' label '_' num2str(re,'%02.f')]);
    mkdir(saveFolder)

    delete(gcp("nocreate"));
    parpool("Processes",8); % initial processing ran on linux computer using 22 cores
    
    parfor e = 1:length(eegFiles)

        % load eeg
        EEG = pop_loadset('filename', eegFiles(e).name,  'filepath', eegFiles(e).folder);
        EEGtmp = EEG;

        % take only the walking and inital standing portions
        EEGtmp = extract_straightSegs(EEGtmp);
        EEGtmp = eeg_checkset(EEGtmp);

        % track inital data duration
        gc_idx1 = checkGC(EEGtmp); % track HS
        starting_gc_num = length(unique(gc_idx1));
        starting_pnts   = EEGtmp.pnts;
        starting_dur    = EEGtmp.pnts/EEGtmp.srate/60;

        % threshold for total amount of data removed: 90% kept
        min_dur = starting_dur*.9;

        %%% remove bad segments using pop_jointprob %%%
        keepLoop = true;
        REJ = 3; % starting jointprob threshold
        EEGtmp2 = eeg_regepochs(EEGtmp,1);
        while keepLoop
            [EEGtmp2_tmp, ~, ~, nrej, rej, com] = pop_jointprobLM(EEGtmp2,1,1:size(EEGtmp2.data,1),REJ,REJ,1,1);
            EEGtmp2_tmp = eeg_epoch2continuous(EEGtmp2_tmp);

            %check
            pnts_removed_jointProb  = EEGtmp.pnts - EEGtmp2_tmp.pnts;
            dur_left = EEGtmp2_tmp.pnts/EEGtmp2_tmp.srate/60;

            if dur_left < min_dur % if more than 90% of gait cycle removed

                if dur_left < (starting_dur*.8)
                    REJ = REJ + 1;
                else
                    REJ = REJ + .5;
                end

                disp('      ---     ');
                disp('Too much data removed! Re-adjusting jointprob threshold...')
                disp(['Total minutes left: ' num2str(dur_left)]);
                disp(['Total minutes left (%): ' num2str(dur_left/starting_dur*100)]);
                disp(['Will try new threshold @ ' num2str(REJ)]);
                disp('      ---     ');
            else
                keepLoop = false;
                epoched_rej = find(rej);
                rm_bounds_jointProb = [epoched_rej*EEG.srate; (epoched_rej*EEG.srate) + (EEG.srate - 1)]';
                sample_mask_jointProb = ones(1,EEGtmp.pnts);
                for x = 1:size(rm_bounds_jointProb,1)
                    sample_mask_jointProb(rm_bounds_jointProb(x,1):rm_bounds_jointProb(x,2)) = 0;
                end
                EEG.etc.jointprob_threshold = REJ;
                EEGtmp2 = [];
                EEGtmp2_tmp = [];
            end
        end

        %%% remove high amp data using clean_windows %%%
        REJ = 5;
        window_crit = 0.3;
        keepLoop = true;
        while keepLoop
            window_tol =  [-Inf REJ];
            [EEGtmp2 sample_mask_cleanWindows]  = clean_windows(EEGtmp,window_crit,window_tol);

            % check removed
            pnts_removed_cleanWindows           = EEGtmp.pnts - EEGtmp2.pnts;
            dur_left = EEGtmp2.pnts/EEGtmp2.srate/60;

            if dur_left < min_dur % if more than 90% of gait cycle removed

                if dur_left < (starting_dur*.8)
                    REJ = REJ + 1;
                else
                    REJ = REJ + .5;
                end

                disp('      ---     ');
                disp('Too much data removed! Re-adjusting clean_windows threshold...')
                disp(['Total minutes left (min): ' num2str(dur_left)]);
                disp(['Total minutes left (%): ' num2str(dur_left/starting_dur*100)]);
                disp(['Will try new threshold @ ' num2str(REJ)]);
                disp('      ---     ');
            else
                keepLoop = false;
                EEG.etc.cleanWindow_threshold = REJ;
                EEGtmp2 = [];
            end
        end

        % merge sample mask and remove bad data seg
        combined_mask = sample_mask_cleanWindows + sample_mask_jointProb;
        combined_mask(combined_mask < 2) = 0;
        combined_mask(combined_mask == 2) = 1;
        keep_intervals = logical2interval(logical(combined_mask));
        EEGtmp = pop_select(EEGtmp, 'point', keep_intervals);

        % store processing info
        EEG.etc.HS_idx_initial          = unique(gc_idx1);
        EEG.etc.HS_idx_passed2ICA       = unique(checkGC(EEGtmp));
        EEG.etc.HS_idx_removed_priorICA = setdiff(EEG.etc.HS_idx_initial , EEG.etc.HS_idx_passed2ICA);
        EEG.etc.DataPercentUsedForICA       = EEGtmp.pnts/starting_pnts;
        EEG.etc.DurationUsedForICA          = EEGtmp.pnts/EEGtmp.srate/60;
        EEG.etc.NumPntsRemoved_cleanWindows        = pnts_removed_cleanWindows;
        EEG.etc.NumPntsRemoved_jointProb           = pnts_removed_jointProb;

        % interpolate, add ref, and re-ref for AMICA
        bad_chans = EEGtmp.chaninfo.removedchans(strcmp({EEGtmp.chaninfo.removedchans.type},'EEG')); % get removed channels
        EEGtmp = pop_interp(EEGtmp, bad_chans, 'spherical'); % interpolate
        EEGtmp.nbchan = EEGtmp.nbchan+1;
        EEGtmp.data(end+1,:) = zeros(1, EEGtmp.pnts);
        EEGtmp.chanlocs(1,EEG.nbchan) = EEGtmp.chaninfo.removedchans(strcmp({EEGtmp.chaninfo.removedchans.type},'REF')); % add back ref channel
        EEGtmp = pop_reref(EEGtmp, []);  % reference to average
        EEGtmp = eeg_checkset(EEGtmp);

        % check rank
        EEGtmp = eeg_checkset(EEGtmp);
        dataRank = getRank(EEGtmp.data);
        %dataRank = sum(eig(cov(double(EEGtmp.data'))) > 1e-7);

        %%% AMICA %%%
        outPath = fullfile(saveFolder,[EEG.setname '_AMICA']);
        weights = []; sphere = []; mods = [];

        % [weights,sphere,mods] = runamica15(EEGtmp.data,'pcakeep',dataRank,...
        %     'max_iter',5000,...
        %     'block_size',256,...
        %     'max_threads',2,...
        %     'fix_init',0,...
        %     'outdir',outPath);

        [weights,sphere,mods] = runamica15(EEGtmp.data,'pcakeep',dataRank,...
            'max_iter',6000,...
            'block_size',256,...
            'max_threads',2,...
            'fix_init',0,...
            'outdir',outPath,...
            'do_reject',1,...
            'numrej',5,...
            'rejstart',2);

        % interpolate, add ref, and re-ref on OG data
        bad_chans = EEG.chaninfo.removedchans(strcmp({EEG.chaninfo.removedchans.type},'EEG')); % get removed channels
        EEG = pop_interp(EEG, bad_chans, 'spherical'); % interpolate
        EEG.nbchan = EEG.nbchan+1;
        EEG.data(end+1,:) = zeros(1, EEG.pnts);
        EEG.chanlocs(1,EEG.nbchan) = EEG.chaninfo.removedchans(strcmp({EEG.chaninfo.removedchans.type},'REF')); % add back ref channel
        EEG = pop_reref(EEG, []);  % reference to average
        EEG = eeg_checkset(EEG);

        % % save interpolated original data to the cleanedDataforICA for relica
    	% EEGtmp.etc = [];
        % EEGtmp.etc.dataRank = dataRank;
        % EEGtmp.etc.original_interpolated_data = EEG.data;
        % EEGtmp.event = [];
        % EEGtmp.urevent = [];
        % tmpFilename = [EEGtmp.subject '_cleanedDataforICA.set'];
        % EEGtmp = pop_saveset(EEGtmp, 'filename', tmpFilename, 'filepath', outPath);

        % add ica weights to OG data
        EEG.icaweights  = weights;
        EEG.icasphere   = sphere;
        EEG.mods        = mods;
        EEG = eeg_checkset(EEG, 'ica');

        % save amica processed data
        EEG = pop_saveset(EEG, 'filename', eegFiles(e).name, 'filepath', saveFolder);
    end
    delete(gcp("nocreate"));
    disp('ICA done.')
end

%% functions
function gc_idx = checkGC(EEG)

event_str = unique({EEG.event.type}); % get gait event labels
% find all events with HS1 in it
if contains(EEG.filepath,'affected','IgnoreCase',true)
    disp('Affected side labels detected, using affected side labels for epoching.')
    HS1_str = event_str(~cellfun('isempty', regexp(event_str, '^MoreAffected_.*_HS.*$')));
    HS1_str(contains(HS1_str,'HS2')) = [];
else
    HS1_str = event_str(contains(event_str,'HS1'));
end

% epoch data
EEG = pop_epoch(EEG,HS1_str, [-1  2], 'epochinfo', 'yes'); % epoch from -1 to 3

gc_idx = NaN(1, length(EEG.epoch));
% check if the epoch has a HS2
for gc = 1:length(EEG.epoch)
    start_idx = find([EEG.epoch(gc).eventlatency{:}] > 0);
    if isempty(start_idx)
        continue
    end
    start_idx = start_idx(1);
    gc_num_cell = EEG.epoch(gc).eventgc_num(start_idx:end);
    gc_num_cell = gc_num_cell(~cellfun(@isempty, gc_num_cell));
    if ~isempty(gc_num_cell)
        gc_idx(gc) = gc_num_cell{1};
    end
end
gc_idx(isnan(gc_idx)) = [];
end