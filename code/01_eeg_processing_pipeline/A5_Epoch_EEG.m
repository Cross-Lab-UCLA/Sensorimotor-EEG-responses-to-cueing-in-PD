%% Gait Cycle Epoching and Time Warping Preparation
%
% This script epochs dataset by gait cycles,
% remove imcomplete cycles and its associated gait data, and extract group
% cycle latencies for future time-warping parameters.
%
% Process all .set files in the selected 04_dipRemoved_*
% folder and saves epoched datasets to 05_epoched_ folder
%
%
% LM 122225
%%
clear all; clc; close all; warning('off','all');
if ispc
    mainDir     = 'C:\Git\DoD-Gait';
else
    mainDir    = '/home/leo/Documents/DoD-Gait';
end
funcpath    = fullfile(mainDir,'code','func');
addpath(funcpath)
disp('Load in 04_dipRemoved folder.')
figDir      = fullfile(mainDir,'reports','data_quality');
workDir     = uigetdir(fullfile(mainDir,'data'),'Load in 04_dipRemoved folder.');

%%
% loop = true;
% while loop
%     % must perform the normal epoching first
%     prompt = "Perform Affected Side Epoching? Y/N [default N]: ";
%     txt = input(prompt,"s");
%     if strcmpi(txt,'Y')
%         PerformAffectedSide = true;
%         loop = false;
%         label       = extractAfter(workDir,'dipRemoved_');
%         saveFolder =  fullfile(mainDir,'data',['05_affectedSideEpoched_' label]);
%         mkdir(saveFolder);
%     elseif strcmpi(txt,'N')
%         PerformAffectedSide = false;
%         loop = false;
%         label       = extractAfter(workDir,'dipRemoved_');
%         saveFolder =  fullfile(mainDir,'data',['05_epoched_' label]);
%         mkdir(saveFolder);
%     end
% end

PerformAffectedSide = false;
label       = extractAfter(workDir,'dipRemoved_');
saveFolder =  fullfile(mainDir,'data',['05_epoched_' label]);
mkdir(saveFolder);
eegFiles    = dir(fullfile(workDir,'*.set'));

%% loop to get time warp values
disp('Calclating time warping values')
gait_events_all = {};
gc_rm_idx = {};
T = [];
removeOutliers = false; % if true, checks each gait cycle or data that is an outlier to the mean gait cycle data

for sub = 1:length(eegFiles) % iteration based on the cluster
    L_TO        = [];
    L_HS        = [];
    R_TO        = [];
    R_HS2       = [];
    no_gc_idx   = [];

    EEG = pop_loadset('filename', eegFiles(sub).name,  'filepath', eegFiles(sub).folder);

    if PerformAffectedSide && contains(EEG.subject,'HC')
        disp('Detected HC when performing affected side epoching. Skipping. ')
        continue
    elseif ~PerformAffectedSide || contains(EEG.subject,'HC')
        disp('normal epoching');
        disp(['Loading ' eegFiles(sub).name])
        event_str = unique({EEG.event.type}); % get gait event labels
        HS1_str = event_str(contains(event_str,'R_HS1'));
        EEG = pop_epoch(EEG,HS1_str, [-1  2], 'epochinfo', 'yes'); % epoch from -1 to 2
        OG_trial_num = EEG.trials;
        [EEG, gait_events_all, gait_table] = getGaitTable(EEG,gait_events_all,sub,removeOutliers);
    else
        disp('Affected side labels detected, using affected side labels for epoching.')
        disp(['Loading ' eegFiles(sub).name])
        event_str = unique({EEG.event.type}); % get gait event labels
        HS1_str = event_str(~cellfun('isempty', regexp(event_str, '^MoreAffected_.*_HS.*$')));
        HS1_str(contains(HS1_str,'HS2')) = [];
        EEG = pop_epoch(EEG,HS1_str, [-1  2], 'epochinfo', 'yes'); % epoch from -1 to 2
        OG_trial_num = EEG.trials;
        [EEG, gait_events_all, gait_table] = getGaitTableAffectedSide(EEG,gait_events_all,sub,removeOutliers);
    end
    
    % check epoch and tables are the same length
    if length(gait_events_all{sub}) ~= EEG.trials
        disp('mismatch between gait cycle used to to calculate warpvals and subject total gait epochs.')
        keyboard
    end
    if length(EEG.epoch) ~= height(gait_table) || ...
            size(EEG.icaact,3) ~= height(gait_table) ||...
            length(gait_events_all{sub}) ~= height(gait_table)
        disp('mismatch between subject gait epochs and gait table.')
        keyboard
    end

    %%% save EEG %%%
    % clean up gait info saved in etc
    EEG.etc.gaitTable = [];
    EEG.etc.turnTable = [];
    EEG.etc.APDMevents = [];
    EEG.etc.APDMraw = [];
    EEG.etc.APDMduration = [];

    % save
    EEG.etc.gait_table = gait_table;
    EEG.etc.gait_events_idx = gait_events_all{sub};
    EEG = eeg_checkset(EEG);
    pop_saveset(EEG, EEG.filename, saveFolder);


    %%% save to table %%%
    T(sub).subject          = EEG.subject;
    T(sub).trialNum_prior   = OG_trial_num;
    T(sub).trialNum_kept    = EEG.trials;
    T(sub).incomplete_trials= EEG.etc.incomplete_cycles_removed;

    if isfield(EEG.etc,'outlier_cycles_removed')
        T(sub).outlier_trials   = EEG.etc.outlier_cycles_removed;
    else
        T(sub).outlier_trials = 0;
    end

end
warpVals = median(cat(1,gait_events_all{:})); % median across all subjects
disp('Epoching is done.');
warning('on','all');
save(fullfile(saveFolder,'warpVals.mat'),"warpVals");

% save table
T = struct2table(T);
writetable(T,fullfile(figDir,['A5_Epoch_' label '.xlsx']));

%%  Functions
function [EEG, gait_events_all, gait_table] = getGaitTable(EEG,gait_events_all,sub,removeOutliers)
no_gc_idx = [];
gait_table = table;
gait_urevent_table = table;
resampled_data_all = nan(EEG.nbchan,200,1);

for m = 1:length(EEG.epoch)
    L_TO_1 = []; L_HS_1 = []; R_TO_1 = []; R_HS2_1 = [];
    L_TO    = cell2mat(EEG.epoch(m).eventlatency(contains(EEG.epoch(m).eventtype, 'L_TO')));
    L_HS    = cell2mat(EEG.epoch(m).eventlatency(contains(EEG.epoch(m).eventtype, 'L_HS')));
    R_TO    = cell2mat(EEG.epoch(m).eventlatency(contains(EEG.epoch(m).eventtype, 'R_TO')));
    R_HS2   = cell2mat(EEG.epoch(m).eventlatency(contains(EEG.epoch(m).eventtype, 'R_HS2')));
    R_HS2_label = EEG.epoch(m).eventtype(contains(EEG.epoch(m).eventtype, 'R_HS2'));
    tmp_table = EEG.epoch(m).eventgait(contains(EEG.epoch(m).eventtype, 'R_HS2'));

    R_HS2_1         = R_HS2(find(R_HS2 > 0, 1, 'first'));
    R_HS2_1_label   = R_HS2_label(find(R_HS2 > 0, 1, 'first'));
    if ~isempty(R_HS2_1)
        L_TO_1    = L_TO(find(L_TO > 0 & L_TO < R_HS2_1, 1, 'first'));
        L_HS_1    = L_HS(find(L_HS > 0 & L_HS < R_HS2_1, 1, 'first'));
        R_TO_1    = R_TO(find(R_TO > 0 & R_TO < R_HS2_1, 1, 'first'));
    end

    if  isempty(L_TO_1) || isempty(L_HS_1) || isempty(R_TO_1) || isempty(R_HS2_1) || ...
            sum(R_TO(R_TO>0)<R_HS2_1) > 1 || R_HS2_1 < R_TO_1 || L_HS_1 < L_TO_1
        % if any of gait cycle label are missing or it is a dt gait cycle
        no_gc_idx = [no_gc_idx; m];
        disp(['Un-usable gait cycle: ' EEG.subject, ' Cycle# = '  num2str(m)])
    elseif  contains(R_HS2_1_label, 'dt')
        no_gc_idx = [no_gc_idx; m];  % store index [ set_number gc_number ]
        disp(['Removing dt cycle: ' EEG.subject, ' Cycle# = '  num2str(m)])
    else
        gait_events_all{sub}(m, :) = [0 L_TO_1 L_HS_1 R_TO_1 R_HS2_1]; % R_HS1 is ALWAYS at zero latency
        T = tmp_table{R_HS2 == R_HS2_1}(1,1:end);
        T.Condition = extractAfter(R_HS2_1_label,'HS2_');
        gait_table         = [gait_table;T];

        % track start and end indexes and latenices from the pre-epoched data
        latencies = cell2mat(EEG.epoch(m).eventlatency);
        R_HS1_event_idx = find(latencies == 0);
        R_HS1_event_idx = R_HS1_event_idx(1);
        R_HS2_event_idx = find(latencies == R_HS2_1);
        R_HS2_event_idx = R_HS2_event_idx(1);
        gait_urevent_table.Condition(m) = T.Condition;
        gait_urevent_table.HS1_idx(m) = EEG.epoch(m).eventurevent{R_HS1_event_idx};
        gait_urevent_table.HS2_idx(m) = EEG.epoch(m).eventurevent{R_HS2_event_idx};
    end

    % get a copy of channel data normalized to the gait cycle
        if ~isempty(R_HS2_1)
            epoch_hs1_idx = find(EEG.times == 0);
            [~, epoch_hs2_idx] = min(abs(EEG.times - R_HS2_1));
            epoch_data = EEG.data(:,epoch_hs1_idx:epoch_hs2_idx,m);
            epoch_times = EEG.times(epoch_hs1_idx:epoch_hs2_idx);
            new_times = linspace(epoch_times(1), epoch_times(end), 200);
            resampled_epoch_data = cell2mat(cellfun(@(row) interp1(epoch_times, row, new_times, 'spline'), ...
                num2cell(epoch_data, 2), 'UniformOutput', false));
            resampled_data_all(:,:,m) = resampled_epoch_data;
        else
            resampled_data_all(:,1:200,m) = nan;
        end

% nexttile()
% plot(epoch_times,epoch_data(1,:));
% nexttile()
% plot(new_times,resampled_epoch_data(1,:))
end

EEG.etc.gaitcycle_StartEnd_idx = gait_urevent_table;

% removed epochs without complete cycles
if ~isempty(no_gc_idx)
    EEG = pop_select(EEG, 'notrial', no_gc_idx);
    EEG = eeg_checkset(EEG);
    gc_rm_idx = find(sum(gait_events_all{sub},2) == 0);
    gait_events_all{sub}(gc_rm_idx,:) = []; % this is needed if a bad trial is at the end;
    resampled_data_all(:,:,no_gc_idx) = [];
    EEG.etc.incomplete_cycles_removed = length(no_gc_idx);
end

% removed epochs that are outliers
if removeOutliers
    rm_outlier_idx = removedOutlierCycles(resampled_data_all);
    if ~isempty(rm_outlier_idx)
        EEG = pop_select(EEG, 'notrial', rm_outlier_idx);
        EEG = eeg_checkset(EEG);
        gait_events_all{sub}(rm_outlier_idx,:) = [];
        gait_table(rm_outlier_idx,:) = [];
        EEG.etc.outlier_cycles_removed = length(rm_outlier_idx);
        EEG.etc.gaitcycle_StartEnd_idx(rm_outlier_idx,:) = [];
    end
end

end

% For affected side epoching
function [EEG, gait_events_all, gait_table] = getGaitTableAffectedSide(EEG,gait_events_all,sub,removeOutliers)
no_gc_idx = [];
gait_table = table;
gait_urevent_table = table;
resampled_data_all = nan(EEG.nbchan,200,1);

for m = 1:length(EEG.epoch)
    L_TO_1 = []; L_HS_1 = []; M_TO_1 = []; M_HS2_1 = [];
    L_TO    = cell2mat(EEG.epoch(m).eventlatency(~cellfun('isempty', regexp(EEG.epoch(m).eventtype, '^LessAffected_.*_TO.*$'))));
    L_HS    = cell2mat(EEG.epoch(m).eventlatency(~cellfun('isempty', regexp(EEG.epoch(m).eventtype, '^LessAffected_.*_HS.*$'))));
    M_TO    = cell2mat(EEG.epoch(m).eventlatency(~cellfun('isempty', regexp(EEG.epoch(m).eventtype, '^MoreAffected_.*_TO.*$'))));
    M_HS2    = cell2mat(EEG.epoch(m).eventlatency(~cellfun('isempty', regexp(EEG.epoch(m).eventtype, '^MoreAffected_.*_HS.*$'))));

    % get the first instance of the event
    M_HS2_1         = M_HS2(find(M_HS2 > 0, 1, 'first'));
    M_HS2_label     = EEG.epoch(m).eventtype(~cellfun('isempty', regexp(EEG.epoch(m).eventtype, '^MoreAffected_.*_HS.*$')));
    M_HS2_1_label   = M_HS2_label{1};

    if ~isempty(M_HS2_1)
        L_TO_1    = L_TO(find(L_TO > 0 & L_TO < M_HS2_1, 1, 'first'));
        L_HS_1    = L_HS(find(L_HS > 0 & L_HS < M_HS2_1, 1, 'first'));
        M_TO_1    = M_TO(find(M_TO > 0 & M_TO < M_HS2_1, 1, 'first'));
    end

    % check to make sure the epoch contains all the stages of a gait cycle
    if  isempty(L_TO_1) || isempty(L_HS_1) || isempty(M_TO_1) || isempty(M_HS2_1) || ...
            sum(M_TO(M_TO>0)<M_HS2_1) > 1 || M_HS2_1 < M_TO_1 || L_HS_1 < L_TO_1
        % if any of gait cycle moment are missing or there are two R TO
        % but only one HS --> this might be from removing noisy epoch
        % that is within a gait cycle
        no_gc_idx = [no_gc_idx; m];  % store index [ set_number gc_number ]
        disp(['Incomplete gait cycle: ' EEG.subject, ' Cycle# = '  num2str(m)])
    elseif  contains(M_HS2_1_label, 'dt')
        no_gc_idx = [no_gc_idx; m];  % store index [ set_number gc_number ]
        disp(['Removing dt cycle: ' EEG.subject, ' Cycle# = '  num2str(m)])
    else
        gait_events_all{sub}(m, :) = [0 L_TO_1 L_HS_1 M_TO_1 M_HS2_1]; % R_HS1 is ALWAYS at zero latency

        % get gait table from the current cycle entire epoch
        gc_idx = find([EEG.epoch(m).eventlatency{:}] >0 & [EEG.epoch(m).eventlatency{:}] <= M_HS2_1);
        tmp_table = EEG.epoch(m).eventgait(gc_idx); % get all gait tables within current epoch
        tmp_table = tmp_table(~cellfun(@isempty, tmp_table)); % remove empty return from indexes without tables
        if isempty(tmp_table) % there shoud be atleast one table per epoch
            keyboard
        end
        tmp_table = tmp_table{end}; % table at the HS2 correspond to the gait data of HS1 to HS2
        tmp_table.Condition = {extractAfter(M_HS2_1_label,'HS1_')};
        if isempty(tmp_table.Condition{1})
            tmp_table.Condition = {extractAfter(M_HS2_1_label,'HS_')};
        end
        gait_table         = [gait_table;tmp_table];

        % track start and end indexes and latenices from the pre-epoched data
        latencies = cell2mat(EEG.epoch(m).eventlatency);
        R_HS1_event_idx = find(latencies == 0);
        R_HS1_event_idx = R_HS1_event_idx(1);
        R_HS2_event_idx = find(latencies == M_HS2_1);
        R_HS2_event_idx = R_HS2_event_idx(1);
        
        gait_urevent_table.Condition(m) = tmp_table.Condition;
        gait_urevent_table.HS1_idx(m) = EEG.epoch(m).eventurevent{R_HS1_event_idx};
        gait_urevent_table.HS2_idx(m) = EEG.epoch(m).eventurevent{R_HS2_event_idx};
    end

    if ~isempty(M_HS2_1)
        epoch_hs1_idx = find(EEG.times == 0);
        [~, epoch_hs2_idx] = min(abs(EEG.times - M_HS2_1));
        epoch_data = EEG.data(:,epoch_hs1_idx:epoch_hs2_idx,m);
        epoch_times = EEG.times(epoch_hs1_idx:epoch_hs2_idx);
        new_times = linspace(epoch_times(1), epoch_times(end), 200);
        resampled_epoch_data = cell2mat(cellfun(@(row) interp1(epoch_times, row, new_times, 'spline'), ...
            num2cell(epoch_data, 2), 'UniformOutput', false));
        resampled_data_all(:,:,m) = resampled_epoch_data;
    else
        resampled_data_all(:,1:200,m) = nan;
    end

end

% removed epochs without complete cycles
if ~isempty(no_gc_idx)
    EEG = pop_select(EEG, 'notrial', no_gc_idx);
    EEG = eeg_checkset(EEG);
    gc_rm_idx = find(sum(gait_events_all{sub},2) == 0);
    gait_events_all{sub}(gc_rm_idx,:) = []; % this is needed if a bad trial is at the end;
    resampled_data_all(:,:,no_gc_idx) = [];
    EEG.etc.incomplete_cycles_removed = length(no_gc_idx);
end
empty_idx = cellfun(@isempty, gait_urevent_table.Condition);
empty_rows = find(empty_idx);
gait_urevent_table(empty_rows,:) = [];
EEG.etc.gaitcycle_StartEnd_idx = gait_urevent_table;

% removed epochs that are outliers
if removeOutliers
    rm_outlier_idx = removedOutlierCycles(resampled_data_all);
    if ~isempty(rm_outlier_idx)
        EEG = pop_select(EEG, 'notrial', rm_outlier_idx);
        EEG = eeg_checkset(EEG);
        gait_events_all{sub}(rm_outlier_idx,:) = [];
        gait_table(rm_outlier_idx,:) = [];
        EEG.etc.outlier_cycles_removed = length(rm_outlier_idx);
        EEG.etc.gaitcycle_StartEnd_idx(rm_outlier_idx,:) = [];
    end
end

end

%% function to removed outlier gait cycles
function rm_idx = removedOutlierCycles(resampled_data_all)

gfp = squeeze(rms(resampled_data_all,1)); % global field power

if any(isnan(gfp))
    keyboard
end

global_threshold1 = median(gfp,2) + std(gfp,0,2)*3; %3 std from mean
global_threshold2 = median(gfp,2) + std(gfp,0,2)*5; %5 std from mean
threshold_percent = 0.05;
std_mask1 = []; std_mask2 = [];

for hs = 1:length(gfp(1,:))
    num_timepoints   = length(gfp(:,1));
    exceed_mask1     = gfp(:,hs) > global_threshold1;
    exceed_count1    = sum(exceed_mask1); 
    std_mask1(hs, :) = exceed_count1 > (threshold_percent * num_timepoints); % Flag if >0.05% exceed

    exceed_mask2     = gfp(:,hs) > global_threshold2;
    exceed_count2    = sum(exceed_mask2); 
    std_mask2(hs, :) = exceed_count2 > 0; 
end

std_mask_combined = std_mask1|std_mask2;
rm_idx = [];
rm_idx = find(std_mask_combined ==1)';

end
