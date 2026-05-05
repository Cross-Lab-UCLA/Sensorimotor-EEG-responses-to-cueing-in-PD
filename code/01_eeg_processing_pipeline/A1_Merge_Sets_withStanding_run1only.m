%% A1_Merge_Sets_withStanding_run1only.m
%
%   This script merges run-1 (or run-2 if run-1 is missing) EEG datasets 
%   for each subject into a single set for clustering and further analysis.
%   Applies filtering, downsampling, and event adjustments for PD subjects 
%   (more/less affected side).
%
%   1. Prompts user to select subjects for analysis.
%   2. Loads and processes walk, walkAuditory, and walkVisual conditions 
%       for each subject. Filter, downsample, remove line noise and 
%       remove high rms chans.
%   3. Merges available conditions into a single EEG set per subject.
%   4. Applies iCanClean and clean_channel
%       - add rng(5489, 'mt19937ar') to clean_channel for parfor processing
%   5. Plots pre/post cleaning results for key channels.
%   6. Updates event labels for affected side in PD subjects.
%   7. Saves merged and cleaned EEG sets and summary tables.
%
% requires plugins iCanClean, clean_rawdata, zapline-plus
%
% LM 091825
%%
clear all; clc; close all;
if ispc
    mainDir     = 'C:\Git\DoD-Gait';
elseif ismac
    mainDir = '/Users/Leo/Git/DoD-Gait';
else
    mainDir = '/home/leo/Documents/DoD-Gait';
end

addpath(fullfile(mainDir,'code','func'))
figDir  = fullfile(mainDir,'reports','data_quality');
mkdir(figDir);
chanDir     = fullfile(figDir,'channelRemovedPlots');
mkdir(chanDir);
plotDir     = fullfile(figDir,'channelCleaningPlots');
mkdir(plotDir);
withStanding = true;

if withStanding
    saveCombDir = fullfile(mainDir,'data','01_combined_withStanding_run1only_moreAffected');
else
    saveCombDir = fullfile(mainDir,'data','01_combined_noStanding_run1only_moreAffected');
end
saveLabel = extractAfter(saveCombDir,'combined_');
dataDir = fullfile(mainDir,'data');
mkdir(saveCombDir);
Ta = readtable(fullfile(dataDir,'affected_side.csv'));

%% add toolboxes
eeglab; close all;

%% subject selection
subjList = dir(fullfile(mainDir,'data','subjects'));
subjList = subjList(cell2mat(cellfun(@(x) contains(x, 'HC', 'IgnoreCase', true) | contains(x, 'PD', 'IgnoreCase', true),...
    {subjList(:).name}, 'UniformOutput', false)));
[idxSubjList,~] = listdlg('PromptString',{'Select Subject(s) for the analysis.',...
    'Single or multiple files can be selected.',''},...
    'SelectionMode','multiple', 'ListSize', [150, 250], 'ListString', {subjList(:).name});

subjList = subjList(idxSubjList);
cellfun(@(x) disp(['>>>>> Selected ' x ' to import <<<<<']), {subjList(:).name}, 'UniformOutput', false);
T1 = cell(length(subjList),1);
T2  = cell2table(cell(length(subjList),12), 'VariableNames', {'subject', 'rest', 'walk', 'walkAuditory', 'walkVisual',...
    'chansRemovedRMS','chansRemovedFlat','chansRemovedCleanChan','chansRemovedTotal','chansKept', 'dataRank', 'totalHS'});
T3 = struct;

%%  start loop to run through each subject 

delete(gcp("nocreate"));
parpool("Processes",7)
parfor s = 1:length(subjList)
    subjID = subjList(s).name(5:end);
    T1_subj_tmp = [];

    T1_subj_tmp.subject = subjList(s).name;
    T2_tmp  = cell(1,12);    % get subject's eeg file names
    T3(s).subject = subjList(s).name;

    subjDir = fullfile(mainDir,'data','subjects',subjList(s).name,'raw');
    eegFiles = dir(fullfile(subjDir,'*.set'));
    T2_tmp{1} = subjID;

    rest = eegFiles(contains({eegFiles.name},'rest','IgnoreCase',true));
    walk1 = eegFiles(contains({eegFiles.name}, 'walk_run-1'));
    walkAudio1 = eegFiles(contains({eegFiles.name}, 'walkAuditory_run-1'));
    walkVisual1 = eegFiles(contains({eegFiles.name}, 'walkVisual_run-1'));
    walk2 = eegFiles(contains({eegFiles.name}, 'walk_run-2'));
    walkAudio2 = eegFiles(contains({eegFiles.name}, 'walkAuditory_run-2'));
    walkVisual2 = eegFiles(contains({eegFiles.name}, 'walkVisual_run-2'));

    %%%%%%%%%%%%% MERGE REST %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    restOn = 0;
    if restOn
        EEG = []; ALLEEG = [];
        for t = 1:length(rest)
            EEG = pop_loadset('filename', rest(t).name,  'filepath', rest(t).folder);
            EEG = cutTrial(EEG,withStanding);
            [EEG, T1_tmp] = cleanTrial(EEG);
            numEvents = length(EEG.event);
            EEG.event(numEvents+1).type = ['START_' EEG.condition] ;
            EEG.event(numEvents+1).latency = 1;
            EEG.event(numEvents+1).urevent = numEvents+1;
            EEG = eeg_checkset(EEG, 'eventconsistency');
            [ALLEEG, ~, ~] = eeg_store(ALLEEG, EEG, t);
        end
        EEG = pop_mergeset(ALLEEG, 1:length(rest), 1);
        EEG.condition = 'restCombined';
        EEG.setname = [EEG.subject '_task-' EEG.condition];
        EEG = pop_saveset(EEG, 'filename', EEG.setname, 'filepath', saveCombDir);
        T2_tmp{2} = length(rest);
    else
        T2_tmp{2} = 'missing';
    end

    %%%%%%%%%%%%% MERGE WALK %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ALLEEG = []; EEG = []; t_counter = 1;
    if ~isempty(walk1)
        disp('processing walk run-1');
        for t = 1:length(walk1)
            EEG = pop_loadset('filename', walk1(t).name,  'filepath', walk1(t).folder);
            EEG = cutTrial(EEG,withStanding);
            [EEG, T1_tmp] = cleanTrial(EEG);
            if ~contains(subjID,'HC')
                EEG = assignAffectedSide(EEG,Ta);
            end
            EEG = assignGaitEvents(EEG);
            [ALLEEG, ~, ~] = eeg_store(ALLEEG, EEG, t_counter);
            t_counter = t_counter + 1;
        end

    elseif ~isempty(walk2) && ~strcmp(subjID,'PD02')
        % PD02 walk_run2 is corrupted - exclude from further analysis

        disp('processing walk run-2');
        for t = 1:length(walk2)
            EEG = pop_loadset('filename', walk2(t).name,  'filepath', walk2(t).folder);
            EEG = cutTrial(EEG,withStanding);
            [EEG, T1_tmp] = cleanTrial(EEG);
            if ~contains(subjID,'HC')
                EEG = assignAffectedSide(EEG,Ta);
            end
            EEG = assignGaitEvents(EEG);
            [ALLEEG, ~, ~] = eeg_store(ALLEEG, EEG, t_counter);
            t_counter = t_counter + 1;
        end
    end

    if isempty(walk1) && isempty(walk2)
        T2_tmp{3} = 'missing';
        T1_subj_tmp.nocue = [];
    else
        T1_subj_tmp.nocue = T1_tmp;

        if ~isempty(walk1)
            T3(s).walk1 = 'x';
        elseif ~isempty(walk2)
            T3(s).walk2 = 'x';
        end
    end

    %%%%%%%%%%%%% MERGE WALK Auditory %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use the first auditory walk trial, unless there is no run-1
    if ~isempty(walkAudio1)
        disp('processing walkAuditory run-1');
        for t = 1:length(walkAudio1)
            EEG = pop_loadset('filename', walkAudio1(t).name,  'filepath', walkAudio1(t).folder);
            EEG = cutTrial(EEG,withStanding);
            [EEG, T1_tmp] = cleanTrial(EEG);
            if ~contains(subjID,'HC')
                EEG = assignAffectedSide(EEG,Ta);
            end
            EEG = assignGaitEvents(EEG);
            [ALLEEG, ~, ~] = eeg_store(ALLEEG, EEG, t_counter);
            t_counter = t_counter + 1;
        end
    elseif ~isempty(walkAudio2)
        disp('processing walkAuditory run-2');
        for t = 1:length(walkAudio2)
            EEG = pop_loadset('filename', walkAudio2(t).name,  'filepath', walkAudio2(t).folder);
            EEG = cutTrial(EEG,withStanding);
            [EEG, T1_tmp] = cleanTrial(EEG);
            if ~contains(subjID,'HC')
                EEG = assignAffectedSide(EEG,Ta);
            end
            EEG = assignGaitEvents(EEG);
            [ALLEEG, ~, ~] = eeg_store(ALLEEG, EEG, t_counter);
            t_counter = t_counter + 1;
        end
    end

    if isempty(walkAudio1) && isempty(walkAudio2)
        T2_tmp{4} = 'missing';
        T1_subj_tmp.audi = [];
    else
        T1_subj_tmp.audi = T1_tmp;
        
        if ~isempty(walkAudio1)
            T3(s).walkAudi1 = 'x';
        elseif ~isempty(walkAudio2)
            T3(s).walkAudi2 = 'x';
        end
    end

    %%%%%%%%%%%%% MERGE WALK Visual %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use the first visual walk trial, unless there is no run-1
    if ~isempty(walkVisual1)
        disp('processing walkVisual run-1');
        for t = 1:length(walkVisual1)
            EEG = pop_loadset('filename', walkVisual1(t).name,  'filepath', walkVisual1(t).folder);
            EEG = cutTrial(EEG,withStanding);
            [EEG, T1_tmp] = cleanTrial(EEG);
            if ~contains(subjID,'HC')
                EEG = assignAffectedSide(EEG,Ta);
            end
            EEG = assignGaitEvents(EEG);
            [ALLEEG, ~, ~] = eeg_store(ALLEEG, EEG, t_counter);
            t_counter = t_counter + 1;
        end
    elseif ~isempty(walkVisual2)
        disp('processing walkVisual run-2');
        for t = 1:length(walkVisual2)
            EEG = pop_loadset('filename', walkVisual2(t).name,  'filepath', walkVisual2(t).folder);
            EEG = cutTrial(EEG,withStanding);
            [EEG, T1_tmp] = cleanTrial(EEG);
            if ~contains(subjID,'HC')
                EEG = assignAffectedSide(EEG,Ta);
            end
            EEG = assignGaitEvents(EEG);
            [ALLEEG, ~, ~] = eeg_store(ALLEEG, EEG, t_counter);
            t_counter = t_counter + 1;
        end
    end

    if isempty(walkVisual1) && isempty(walkVisual2)
        T2_tmp{5} = 'missing';
        T1_subj_tmp.visual = [];
    else
        T1_subj_tmp.visual = T1_tmp;
        if ~isempty(walkVisual1)
            T3(s).walkVisual1 = 'x';
        elseif ~isempty(walkVisual2)
            T3(s).walkVisual2 = 'x';
        end
    end

    %%% Merge all condition sets  %%%
    for tt = 1:length(ALLEEG)-1
        if tt == 1
            EEG_A = ALLEEG(tt);
        else
            EEG_A = EEG_merged; % the previous merged EEG
        end
        EEG_B = ALLEEG(tt+1);

        labelsA = {EEG_A.chanlocs.labels};
        labelsB = {EEG_B.chanlocs.labels};
        missing_in_A = setdiff(labelsB, labelsA); %labels in B but not in A:
        missing_in_B = setdiff(labelsA, labelsB);

        if ~isempty(missing_in_B)
            EEG_A = pop_select(EEG_A,'nochannel', missing_in_B);
        end
        if ~isempty(missing_in_A)
            EEG_B = pop_select(EEG_B,'nochannel', missing_in_A);
        end
        EEG_merged = pop_mergeset(EEG_A, EEG_B, 1);
    end
    EEG = EEG_merged;
    EEG_merged = []; ALLEEG = []; EEG_A = []; EEG_B = [];
    
    %%% plot pre ICC %%%
    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf, 'Color', 'w');
    COI = {'3LA' '3RA'};
    COI_idx = find(ismember({EEG.chanlocs.labels},COI));
    tl = tiledlayout(2,length(COI_idx));
    title(tl,'Pre vs Post ICCleaning');
    for c = 1:length(COI_idx)
        nexttile
        h1 = plot(EEG.times/1000/60,EEG.data(COI_idx(c),:),'LineWidth',1.5); hold on
        title([COI{c} ' - timeseries'])
        xlabel('Time (minute)'); ylabel('Amplitude (µV)')
        grid on

        nexttile
        [pxx,f] = pwelch(EEG.data(COI_idx(c),:),EEG.srate,EEG.srate/2,EEG.srate,EEG.srate);
        h2 = plot(f(1:70),10*log10(pxx(1:70))); hold on
        title([COI{c} ' - PSD'])
        xlabel('Frequency (Hz)'); ylabel('Power (dB)')
        grid on
    end

    %%% ICC %%%
    % apply iCanClean parameters
    params = [];
    params.rhoSqThres_source = .8;
    params.filtYtype = 'Notch';
    params.filtYfreq = [4 50];
    params.plotStatsOn = false;
    %params.rerefX = 'yes-fullrank';
    if ~contains(EEG.condition,'rest')
        EEG = iCanClean(EEG, [1:length(EEG.chanlocs)], [1:length(EEG.chanlocs)] ,0, params);
    end

    %%% remove channels using clean_channels %%%
    rms_removed_chans = {EEG.chaninfo.removedchans(strcmp({EEG.chaninfo.removedchans.type},'EEG')).labels};
    [EEG,~] = clean_channels(EEG,.75,[],[],[],150,[]); % corr_threshold @ .75; mannually set rng for clean_channels at line 114 to rng(5489, 'mt19937ar')
    bad_chans = EEG.chaninfo.removedchans(strcmp({EEG.chaninfo.removedchans.type},'EEG'));
    if ~isempty(bad_chans)
        cleanraw_removed_ch = {bad_chans.labels};
        cleanraw_removed_ch = setdiff(cleanraw_removed_ch, rms_removed_chans, 'stable');
    else
        cleanraw_removed_ch = [];
    end

    % plot removed channels
    fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
    topoplot([],EEG.chanlocs, 'style', 'blank', 'drawaxis', 'on', 'electrodes', ...
        'labelpoint', 'chaninfo', EEG.chaninfo); hold on
    saveFigName = fullfile(chanDir,[EEG.subject '.png']);
    saveas(fig2,saveFigName);
    close(fig2)
 
    % update name
    EEG.condition = 'walkCombined';
    EEG.setname = [EEG.subject '_task-' EEG.condition];

    %%% interpolate channels %%%
    % bad_chans = EEG.chaninfo.removedchans(strcmp({EEG.chaninfo.removedchans.type},'EEG'));
    % EEG = pop_interp(EEG, bad_chans, 'spherical'); 
    %%% add back ref channel %%%
    % EEG.nbchan = EEG.nbchan+1;
    % EEG.data(end+1,:) = zeros(1, EEG.pnts);
    % EEG.chanlocs(1,EEG.nbchan) = EEG.chaninfo.removedchans(strcmp({EEG.chaninfo.removedchans.type},'REF'));
    %%% re-reference %%%
    % EEG = pop_reref(EEG, []);  % reference to average
    % EEG = eeg_checkset(EEG);
    
    % check rank
    EEG.etc.dataRank = getRank(EEG.data);

    % assign gait cycle numbers
    HS_idx = find(contains({EEG.event.type},'R_HS2'));
    counter = 0;
    for hs = 1:length(HS_idx)
        counter = counter + 1;
        EEG.event(HS_idx(hs)).gc_num = counter;
    end

    % track number of chans kept
    EEG.etc.bad_chans_A1 = bad_chans;
    fn = fieldnames(T1_subj_tmp); % track down flatline chans
    flat_removed_ch = {};
    for k = 1:numel(fn)
        thisField = T1_subj_tmp.(fn{k});
        if isstruct(thisField) && isfield(thisField,'channels_removedbyFlatline')
            flat_removed_ch = [flat_removed_ch, thisField.channels_removedbyFlatline];
        end
    end
    flat_removed_ch = unique(flat_removed_ch);

    T2_tmp{6} = length(rms_removed_chans) - length(flat_removed_ch);
    T2_tmp{7} = length(flat_removed_ch);
    T2_tmp{8} = length(cleanraw_removed_ch);
    T2_tmp{9} = length(bad_chans);
    T2_tmp{10} = EEG.nbchan;
    T2_tmp{11} = EEG.etc.dataRank;
    T2_tmp{12} = length(HS_idx);
    T2{s,:} = T2_tmp;
    
    T1_subj_tmp.channelRemovedByRMS         = T2_tmp{6};
    T1_subj_tmp.channelRemovedByFlatline    = T2_tmp{7};
    T1_subj_tmp.channelRemovedByCleanChan   = T2_tmp{8};
    T1_subj_tmp.channelRemovedTotal         = T2_tmp{9};
    T1_subj_tmp.channelKept                 = T2_tmp{10};
    T1_subj_tmp.dataRank                    = T2_tmp{11};
    T1_subj_tmp.totalHS                     = T2_tmp{12};

    %%% plot post ICC %%%
    COI = {'3LA' '3RA'};
    COI_idx = find(ismember({EEG.chanlocs.labels},COI));
    ax = findall(f1,'type','axes'); % get all current axes
    ax = flipud(ax); % tiledlayout returns reversed order
    for c = 1:length(COI_idx)
        axes(ax((c-1)*2+1))
        plot(EEG.times/1000/60,EEG.data(COI_idx(c),:),'LineWidth',1);hold on
        title([COI{c} ' - timeseries']);
        legend({'Pre-clean','Post-clean'},'Location','best')

        axes(ax((c-1)*2+2))
        [pxx,f] = pwelch(EEG.data(COI_idx(c),:),EEG.srate,EEG.srate/2,EEG.srate,EEG.srate);hold on
        plot(f(1:70),10*log10(pxx(1:70)))
        title([COI{c} ' - PSD'])
        legend({'Pre-clean','Post-clean'},'Location','best')
    end
    saveFigName = fullfile(plotDir,[EEG.subject '.png']);
    saveas(f1,saveFigName);
    close(f1)

    % save
    EEG = eeg_checkset(EEG);
    EEG = pop_saveset(EEG, 'filename', EEG.setname, 'filepath', saveCombDir);
    T1{s} = T1_subj_tmp;
end
delete(gcp("nocreate"));

% write to reports\data_quality
writetable(T2,fullfile(figDir,['A1_' saveLabel '.xlsx']));
filename = fullfile(figDir, ['A1_trackingStruct_' saveLabel '.mat']);
save(filename,'T1')
writetable(struct2table(T3),fullfile(figDir,['A1_TrialsUsed.xlsx']));
disp('done');


%%
%% FUNCTIONS
%%
function EEG = cutTrial(EEG,withStanding)
% ID and remove data 2 before 1st heel strike and 2 second after last heel strike
if contains(EEG.condition,'rest')
    startLatency = EEG.event(contains({EEG.event.type},'apdm_start')).latency;
    endLatency = EEG.event(contains({EEG.event.type},'apdm_end')).latency;
    EEG = pop_select(EEG,'point', [startLatency endLatency]);
else
    if withStanding
        HS1_lat         = [EEG.event(find(contains({EEG.event(:).type}, 'standing_start'))).latency];
    else
        HS1_lat         = [EEG.event(find(contains({EEG.event(:).type}, 'HS1'))).latency];
    end
    HS2_lat         = [EEG.event(find(contains({EEG.event(:).type}, 'HS2'))).latency];
    startLatency    = HS1_lat(1);
    endLatency      = HS2_lat(end) + 2*EEG.srate;     % last HS plus 2 seconds after it
    EEG = pop_select(EEG,'point', [startLatency endLatency]);
end
end

%%
function EEG = assignGaitEvents(EEG)

hs_idx = find(contains({EEG.event.type},{'HS' 'standing_start' 'standing_end'}));
for h = 1:length(hs_idx)
    EEG.event(hs_idx(h)).type = [EEG.event(hs_idx(h)).type '_' EEG.condition]; % add condition label to events
end
numEvents = length(EEG.event);
EEG.event(numEvents+1).type = ['START_' EEG.condition] ;
EEG.event(numEvents+1).latency = 1;
EEG.event(numEvents+1).urevent = numEvents+1;

for i = 1:numel(EEG.event)
    EEG.event(i).condition = EEG.condition; %
end

EEG = eeg_checkset(EEG, 'eventconsistency');

end

%%
function EEG = assignAffectedSide(EEG,Ta)
% change HS label from R and L to affected

for ev = 1:length(EEG.event)
    if contains(EEG.event(ev).type,'R_HS2')
        % find the index of the L_HS that came immediately before the
        % current R_HS2
        prev_L_HS_idx = find(contains({EEG.event(1:ev-1).type}, 'L_HS'), 1, 'last');

        % make sure there is no extra HS2 between the L_HS and R_HS2
        if ~any(contains({EEG.event(prev_L_HS_idx:ev-1).type},'R_HS2'))

            % assign gait table to L heelstrike events
            EEG.event(prev_L_HS_idx).gait = EEG.event(ev).gait;
            originalLabel = EEG.event(prev_L_HS_idx).type;
            EEG.event(prev_L_HS_idx).type = strrep(originalLabel, 'HS_dt', 'HS1_dt');
        end
    end
end

% reassigned HS labels to more affected side
sub_idx = find(strcmp(Ta.subject,EEG.subject));
if Ta.affected_side_label{sub_idx} == 'R'
    for i = 1:numel(EEG.event)
        if contains(EEG.event(i).type, 'R_')
            EEG.event(i).type = strrep(EEG.event(i).type, 'R_', 'MoreAffected_R_');
        elseif contains(EEG.event(i).type, 'L_')
            EEG.event(i).type = strrep(EEG.event(i).type, 'L_', 'LessAffected_L_');
        end
    end
elseif Ta.affected_side_label{sub_idx} == 'L'
    for i = 1:numel(EEG.event)
        if contains(EEG.event(i).type, 'L_')
            EEG.event(i).type = strrep(EEG.event(i).type, 'L_', 'MoreAffected_L_');
        elseif contains(EEG.event(i).type, 'R_')
            EEG.event(i).type = strrep(EEG.event(i).type, 'R_', 'LessAffected_R_');
        end
    end
end

end