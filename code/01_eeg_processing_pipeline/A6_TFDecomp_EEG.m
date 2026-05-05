%% Time-Frequency Decomposition for Gait Cycle Analysis
%
% This script performs time-frequency decomposition for ICA components 
% from epoched EEG data. It uses EEGLAB's `newtimef` with 
% time-warping to align gait cycles across all subjects.
%
% Prompts user to select a folder containing epoched EEG datasets (05_epoched_*).
% Loads precomputed gait cycle time-warping values (`warpVals.mat`) from
% previous step. Remove outlier epochs using modified jointprob function.
% Runs `newtimef` to compute time-frequency data for each ICA component in each subject.
% Applies time-warping to align gait events across trials (each gait cycle).
% Normalize ERSP using a custom function normalizeERSP.
% Saves ERSP data for each subject to a '06_ERSP_*' folder.
% 
% LM 122225
%%
clear all; 
if ispc
    mainDir     = 'C:\Git\DoD-Gait';
elseif ismac
    mainDir = '/Users/Leo/Git/DoD-Gait';
else
    mainDir    = '/home/leo/Documents/DoD-Gait';
end

% load in libraries
funcpath   = fullfile(mainDir,'code','func');
addpath(funcpath)
eeglab;
ft_defaults;
clc; close all

% set up folder
dataDir = fullfile(mainDir,'data');
figDir      = fullfile(mainDir,'reports','data_quality');

% prompt user for folder
disp('Load in 05_epoched folder.')
path = uigetdir(dataDir,'Load in 05_epoched folder.'); % select folder containing epoched dataset
label = extractAfter(path,'epoched_');
saveSubjFolder = fullfile(dataDir,['06_ERSP_' label]); % set up folder for saving time-frequency data
mkdir(saveSubjFolder)
load(fullfile(path,'warpVals.mat'));

% get files
files = dir(fullfile(path,'*.set'));
%% run newtimef
% open parallel processing
cy = [2 0.5];
EEG = [];
delete(gcp("nocreate"));
parpool("Processes",7)

% run loop for ersp extraction
parfor s = 1:length(files)
    fprintf('Processing subject %s\n', files(s).name);
    filename = fullfile(files(s).folder,files(s).name);
    EEG = load(filename,"-mat");
    EEG = eeg_checkset(EEG,'ica');

    % save to ERSP
    ERSP            = [];
    ERSP.subject    = EEG.subject;
    ERSP.group      = EEG.group;
    ERSP.dipfit     = EEG.dipfit;

    % track gc indexes this stage
    gc_kept = [];
    for gc = 1:length(EEG.epoch)
        start_idx = find([EEG.epoch(gc).eventlatency{:}] > 0);
        start_idx = start_idx(1);
        gc_num_cell = EEG.epoch(gc).eventgc_num(start_idx:end);
        gc_num_cell = gc_num_cell(~cellfun(@isempty, gc_num_cell));
        gc_kept(gc) = gc_num_cell{1};
    end
    missingGC = setdiff(EEG.etc.HS_idx_initial, gc_kept);
    ERSP.HS_idx.inital = EEG.etc.HS_idx_initial;
    ERSP.HS_idx.removed4ICA = EEG.etc.HS_idx_removed_priorICA;
    ERSP.HS_idx.removedAfterEpoch = missingGC;
    ERSP.HS_idx.keptAfterEpoch = gc_kept;

    % perform jointprob to remove outliers
    REJ = 3; % threshold for jointprob
    % check outliers for both channels and ICs
    [~, ~, ~, ~, removedIdx1] = pop_jointprobLM(EEG,1,1:size(EEG.data,1),REJ,REJ,1,1);
    [~, ~, ~, ~, removedIdx2] = pop_jointprobLM(EEG,0,1:size(EEG.icaact,1),REJ,REJ,1,1);
    
    removedIdx = unique([find(removedIdx1) find(removedIdx2)]);

    ERSP.HS_idx.removed4TF  = gc_kept(removedIdx);
    ERSP.HS_idx.kept4TF     = gc_kept;
    ERSP.HS_idx.kept4TF(removedIdx) = [];

    % removed bad epochs 
    EEG = pop_select(EEG, 'notrial', removedIdx);
    EEG = eeg_checkset(EEG);
    EEG.etc.gait_events_idx(removedIdx,:) = [];
    EEG.etc.gait_table(removedIdx,:) =[];
    
    % save some RAM
    EEG.data = [];
    EEG.event = [];
    EEG.urevent = [];

    % run newtimef
    base = []; tfdata_all = []; new_times   = [];
    for ic = 1:size(EEG.icaact,1)
        %running newtimef to get the baseline (powbase) from combined conditions
        [~, ~, ~, times , freqs , ~, ~, tfdata] = ...
            newtimef(EEG.icaact(ic,:,:), EEG.pnts,...
            [EEG.times(1) EEG.times(end)], EEG.srate, ...
            'cycles', cy,...
            'freqs', [4 50],...
            'padratio', 2, ...
            'timesout',-1,...
            'timewarp', EEG.etc.gait_events_idx, 'timewarpms', warpVals,...
            'baseline',[warpVals(1) warpVals(end)],...
            'plotersp','off', ...
            'plotitc', 'off');

        HS1 = find(times == 0);
        HS2 = find(times >= warpVals(end));
        HS2 = HS2(1);
        new_times = times(HS1:HS2);
        % data = bsxfun(@rdivide, tf, base);
        tfdata_all(ic,:,:,:) = tfdata(:,HS1:HS2,:);
    end

    % get baseline from all conditions
    data = tfdata_all.*conj(tfdata_all);
    mean_data = mean(data,4);
    base = mean(mean_data,3);

    % save tf data to ERSP struct
    ERSP.times           = new_times;
    ERSP.freqs           = freqs;
    ERSP.data            = tfdata_all;
    ERSP.baselineAllCond = base;
    ERSP.warpVals        = warpVals;
    ERSP.gait_event_idx  = EEG.etc.gait_events_idx;
    ERSP.gait_table      = EEG.etc.gait_table;
    ERSP.chanlocs        = EEG.chanlocs;
    ERSP.icawinv         = EEG.icawinv;

    % clear for RAM
    EEG = []; tfdata_all = []; data = []; mean_data = [];

    % % remove outlier
    % ERSP = removeOutliersERSP(ERSP,removalFigDir);
    % for tt = 1:length(ERSP.keep)
    %     ERSP.HS_idx.removedAfterOutlier{tt} = ERSP.HS_idx.keptAfterEpoch(~ERSP.keep{tt});
    %     ERSP.HS_idx.keptAfterOutlier{tt} = ERSP.HS_idx.keptAfterEpoch(ERSP.keep{tt});
    % end

    % normalize ERSP
    ERSP.tfdata = normalizeERSP(ERSP);
    ERSP = rmfield(ERSP, 'data'); % clear for space

    % save EEG struct to file
    filename = [ERSP.subject '_ersp_gait.mat'];
    fprintf('Saving subject %s\n', ERSP.subject);    
    parsaveERSP(fullfile(saveSubjFolder,filename),ERSP)
    fprintf('Finished subject %s\n', ERSP.subject);
end

%close parallel processing
delete(gcp('nocreate'))
disp('newtimef done.');

function parsaveERSP(filename,ERSP)
    save(filename, '-struct', 'ERSP','-v7.3','-nocompression');
end