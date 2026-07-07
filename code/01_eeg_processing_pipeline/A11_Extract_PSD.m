%% A11_Extract_PSD
%   Extracts power spectral density for the one selected IC per subject per ROI,
%   separately for walking and standing, and decomposes each spectrum with FOOOF
%   into aperiodic and periodic components. Can be run after A7.
%
%   Requires python set up for matlab (the fooof_mat wrapper calls it).
%
%   Labels:
%   spectraAll_*           retained gait-cycle segments only, pooled across the
%                          three walking conditions (1 row). Excludes standing
%                          and any data not kept for time-frequency analysis.
%   spectra_*              walking, per condition (3 rows)
%   spectra_standing_*     standing, per condition (3 rows)
%
%   output labels
%   _db                    10*log10(power), dB - as returned by spectopo
%   _power                 linear power, uV^2/Hz - from 10.^(dB/10)
%   _percent               percent of total power across 4-50 Hz
%   _fooof_log10           FOOOF power_spectrum - ap_fit = log10(total / aperiodic),
%                          a dimensionless log ratio in bels. This is the flattened
%                          spectrum (peaks plus model error), not the fitted peaks.
%   _fooof_db              the same log ratio x10, in dB
%   _fooof_power           total - aperiodic in linear power, uV^2/Hz. A difference,
%                          not a unit variant of _fooof_log10, and negative wherever
%                          the aperiodic fit sits above the data.
%   ap_fit_log10           FOOOF aperiodic fit, log10 power
%   ap_fit_db              the same quantity x10, in dB
%   ap_fit_power           the same quantity as 10.^ap_fit_log10, uV^2/Hz
%
%
% LM 091526
%%
clear all; clc; eeglab; close all
if ispc
    mainDir     = 'C:\clab\DoD-Gait';
else
    mainDir    = '/home/leo/Documents/DoD-Gait';
end
funcPath    = fullfile(mainDir,'code','01_eeg_processing_pipeline','func');
addpath(funcPath);
addpath(fullfile(funcPath,'fooof_mat'))
dataDir = fullfile(mainDir,'data');
figDir      = fullfile(mainDir,'reports','data_quality');

% prompt user for eeg folder
disp('Load in 04_dipRemoved folder.')
workDir     = uigetdir(fullfile(mainDir,'data'),'Load in 04_dipRemoved* folder.'); % user select folder (i.e., 03_dipfitted)
%label       = extractAfter(workDir,'dipRemoved_');
eegFiles    = dir(fullfile(workDir,'*.set'));

% create save dir
saveDir = fullfile(dataDir,'08_PSD');
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

%% loop to get time PSD
table1 = table(); % table for tracking
counter = 0;
PSD = [];
f_range = [4 50];

for sub = 1:length(eegFiles) % iteration based on the cluster

    % load in eeg file
    EEG = pop_loadset('filename', eegFiles(sub).name,  'filepath', eegFiles(sub).folder);
    
    % use selected subjects and dipoles only
    erspFolder = fullfile(dataDir,'07_selected');
    erspFiles = dir(fullfile(erspFolder,['*' EEG.subject '*']));
    
    for es = 1:length(erspFiles)
        ic = extractBefore(erspFiles(es).name,'_sub');
        load(fullfile(erspFiles(es).folder,erspFiles(es).name))

        % match the location of dipole to track idx
        dipole_pos = tfdata.dipfit.posxyz;
        all_pos = vertcat(EEG.dipfit.model.posxyz);
        diffs = abs(all_pos - dipole_pos);
        match_idx = find(sum(diffs,2) < 0.001);
        if isempty(match_idx) || length(match_idx) > 1
            keyboard
        end

        % split standing
        standing_start  = find(contains({EEG.event.type}, 'standing_start'));
        standing_start_lat = [EEG.event(standing_start).latency];
        standing_end    = find(contains({EEG.event.type}, 'standing_end'));
        standing_end_lat = [EEG.event(standing_end).latency];

        trial_types = {EEG.event(standing_start).type};
        expected_trials = {'walk-run', 'walkAuditory-run', 'walkVisual-run'};

        spectra_standing = nan(3,257);
        for ss = 1:length(standing_start_lat)
            % slot by condition name so standing rows align with the walking
            % rows below (1=walk, 2=auditory, 3=visual), independent of order
            row = find(cellfun(@(t) contains(trial_types{ss}, t), expected_trials));
            if isempty(row)
                keyboard
            end
            EEGtmp = pop_select(EEG, 'point', [standing_start_lat(ss) standing_end_lat(ss)]);
            EEGtmp = eeg_checkset(EEGtmp);
            data_standing = EEGtmp.icaact(match_idx,:);
            [spectra_standing(row,:),freqs] = spectopo(data_standing,0, EEGtmp.srate,...
                'winsize', EEGtmp.srate, 'nfft', EEGtmp.srate*2, ...
                'overlap',EEGtmp.srate*.5, 'plot','off');
        end

        % split EEG by conditions
        start_events = find(contains({EEG.event.type}, 'START_'));
        start_events_lat = [EEG.event(start_events).latency];
        end_events_lat = [start_events_lat(2:end)-1, EEG.pnts];

        keepIdx = tfdata.HS_idx.kept4TF;
        keepArray = nan(length(keepIdx),2);

        for hs = 1:length(keepIdx)
            curr_hs_idx = find(cellfun(@(x) isequal(x, keepIdx(hs)), {EEG.event.gc_num}));
            RTO_idx = find(contains({EEG.event(1:curr_hs_idx).type},'R_TO'));
            HS1_idx = find(contains({EEG.event(1:RTO_idx(end)).type},'R_HS1'));
            HS1_lat = EEG.event(HS1_idx(end)).latency;
            HS2_lat = EEG.event(curr_hs_idx).latency;
            if HS2_lat - HS1_lat > 2000
                keyboard
            end
            keepArray(hs,:) = [HS1_lat HS2_lat];
        end

        % merge times that are consective
        collapsed = keepArray(1, :);
        for i = 2:size(keepArray, 1)
            % If current start equals the previous end, extend the interval
            if abs(keepArray(i,1) - collapsed(end,2)) < 1e-10
                collapsed(end,2) = keepArray(i,2);  % Extend the interval
            else
                collapsed = [collapsed; keepArray(i,:)];  % Start a new one
            end
        end
        keepArrayCollapsed = [floor(collapsed(:,1)) ceil(collapsed(:,2))];
        
        % create EEG with kept HS arrays
        EEGtmp = pop_select(EEG, 'point', keepArrayCollapsed);    
        EEGtmp = eeg_checkset(EEGtmp);

        % set up subject psd matrix
        spectra = nan(3,257);
        
        % get no-cue PSD
        nocue1 = find(contains({EEGtmp.event.type}, 'R_HS1_walk-run'));
        nocue1 = floor(EEGtmp.event(nocue1(1)).latency);
        nocue2 = find(contains({EEGtmp.event.type}, 'R_HS2_walk-run'));
        nocue2 = ceil(EEGtmp.event(nocue2(end)).latency);
        EEG_nocue = pop_select(EEGtmp, 'point', [nocue1 nocue2]);
        data_nocue = EEG_nocue.icaact(match_idx,:);
        boundary_events = strcmp({EEG_nocue.event.type}, 'boundary');
        boundary_samples = round([EEG_nocue.event(boundary_events).latency]);
        boundary_samples(end) = []; % last boundary is after the last data point
        [spectra(1,:),freqs] = spectopo(data_nocue,0, EEG_nocue.srate,...
            'winsize', EEG_nocue.srate, 'nfft', EEG_nocue.srate*2,...
            'overlap',EEGtmp.srate*.5,...
            'boundaries',boundary_samples,'plot','off');

        % get audi PSD
        audi1 = find(contains({EEGtmp.event.type}, 'R_HS1_walkAuditory'));
        audi1 = floor(EEGtmp.event(audi1(1)).latency);
        audi2 = find(contains({EEGtmp.event.type}, 'R_HS2_walkAuditory'));
        audi2 = ceil(EEGtmp.event(audi2(end)).latency);
        EEG_audi = pop_select(EEGtmp, 'point', [audi1 audi2]);
        data_audi = EEG_audi.icaact(match_idx,:);
        boundary_events = strcmp({EEG_audi.event.type}, 'boundary');
        boundary_samples = round([EEG_audi.event(boundary_events).latency]);
        boundary_samples(end) = []; % last boundary is after the last data point
        [spectra(2,:),freqs] = spectopo(data_audi,0, EEG_audi.srate,...
            'winsize', EEG_audi.srate, 'nfft', EEG_audi.srate*2,...
            'overlap',EEGtmp.srate*.5,...
            'boundaries',boundary_samples,'plot','off');

        % get audi PSD
        visual1 = find(contains({EEGtmp.event.type}, 'R_HS1_walkVisual-run'));
        if ~isempty(visual1)
        visual1 = floor(EEGtmp.event(visual1(1)).latency);
        visual2 = find(contains({EEGtmp.event.type}, 'R_HS2_walkVisual-run'));
        visual2 = ceil(EEGtmp.event(visual2(end)).latency);
        EEG_visual = pop_select(EEGtmp, 'point', [visual1 visual2]);
        data_visual = EEG_visual.icaact(match_idx,:);
        boundary_events = strcmp({EEG_visual.event.type}, 'boundary');
        boundary_samples = round([EEG_visual.event(boundary_events).latency]);
        boundary_samples(end) = []; % last boundary is after the last data point
        [spectra(3,:),freqs] = spectopo(data_visual,0, EEG_visual.srate,...
            'winsize', EEG_visual.srate, 'nfft', EEG_visual.srate*2,...
            'overlap',EEGtmp.srate*.5,...
            'boundaries',boundary_samples,'plot','off');
        end

        % get all PSD
        data_all = EEGtmp.icaact(match_idx,:);
        boundary_events = strcmp({EEGtmp.event.type}, 'boundary');
        boundary_samples = round([EEGtmp.event(boundary_events).latency]);
        boundary_samples(end) = []; % last boundary is after the last data point
        [spectraAll,freqs] = spectopo(data_all,0, EEGtmp.srate,...
            'winsize', EEGtmp.srate, 'nfft', EEGtmp.srate*2,...
            'overlap',EEGtmp.srate*.5,...
            'boundaries',boundary_samples,'plot','off');

        counter = counter + 1;
        PSD(counter).subject    = tfdata.subject;
        PSD(counter).roi        = tfdata.roi;
        PSD(counter).group      = tfdata.group;
        PSD(counter).subgroup   = tfdata.subgroup{1};
        PSD(counter).age        = tfdata.age;
        PSD(counter).height     = tfdata.height;
        PSD(counter).weight     = tfdata.weight;
        PSD(counter).moca       = tfdata.MOCA;
        PSD(counter).updrs3     = tfdata.UPDRS3;
        PSD(counter).dipfit     = tfdata.dipfit;
        PSD(counter).labels     = tfdata.labels;
        PSD(counter).gait_table = tfdata.gait_table;
        PSD(counter).freqs      = freqs;
        
        PSD(counter).spectraAll_db = spectraAll;
        PSD(counter).spectra_db    = spectra;

        % convert to power uV^2/Hz
        spectraAll_power            = 10.^(spectraAll/10);
        spectra_power               = 10.^(spectra/10);
        spectra_standing_power      = 10.^(spectra_standing/10);

        PSD(counter).spectra_power          = spectra_power;
        
        % normalize by total power
        f1 = find(freqs == f_range(1));
        f2 = find(freqs == f_range(2));
        PSD(counter).freqs_percent = freqs(f1:f2);
        PSD(counter).spectraAll_percent = spectraAll_power(f1:f2)./sum(spectraAll_power(f1:f2)) * 100;
        PSD(counter).spectra_percent = spectra_power(:,f1:f2)./sum(spectra_power(:,f1:f2),2) * 100;

        % fooof
        settings = struct(); % Default FOOOF settings
        settings.peak_width_limits = [2 12]; % set the minimum band-width at twice the freq res (i.e., if bin width is 1Hz, set min to 2)
        fooof_data = fooof(freqs, spectraAll_power, f_range, settings, true);
        PSD(counter).spectraAll_fooof_log10 = fooof_data.power_spectrum - fooof_data.ap_fit;
        ap_fit_all = fooof_data.ap_fit;
        spectra_fooof   = nan(3,length(PSD(counter).spectraAll_fooof_log10));
        ap_fit          = nan(3,length(PSD(counter).spectraAll_fooof_log10));
        for ff = 1:size(spectra_power,1)
            if any(isnan(spectra_power(ff,:)))
                continue
            end
            fooof_data          = fooof(freqs, spectra_power(ff,:), f_range, settings, true);
            spectra_fooof(ff,:) = fooof_data.power_spectrum - fooof_data.ap_fit;
            ap_fit(ff,:)        = fooof_data.ap_fit;
        end
        
        PSD(counter).freqs_fooof    = fooof_data.freqs;
        PSD(counter).spectra_fooof_log10 = spectra_fooof; % OG spectra - ap fit
        PSD(counter).ap_fit_log10        = ap_fit;

        % dB versions; fooof works in log10 power, so x10 converts to dB
        PSD(counter).spectraAll_fooof_db = 10*PSD(counter).spectraAll_fooof_log10;
        PSD(counter).spectra_fooof_db    = 10*spectra_fooof;
        PSD(counter).ap_fit_db           = 10*ap_fit;

        % linear power versions (uV^2/Hz); the periodic part must be subtracted in
        % linear space, since 10.^spectra_fooof_log10 is a ratio rather than power
        PSD(counter).ap_fit_power           = 10.^ap_fit;
        PSD(counter).spectra_fooof_power    = spectra_power(:,f1:f2)  - 10.^ap_fit;
        PSD(counter).spectraAll_fooof_power = spectraAll_power(f1:f2) - 10.^ap_fit_all;

        % standing PSD
        PSD(counter).spectra_standing_db = spectra_standing;
        PSD(counter).spectra_standing_power = spectra_standing_power;
        PSD(counter).spectra_standing_percent = spectra_standing_power(:,f1:f2)./sum(spectra_standing_power(:,f1:f2),2) * 100;
        spectra_standing_fooof  = nan(3,length(PSD(counter).spectraAll_fooof_log10));
        ap_fit_standing         = nan(3,length(PSD(counter).spectraAll_fooof_log10));

        for ff = 1:size(spectra_standing_power,1)
            if any(isnan(spectra_standing_power(ff,:)))
                continue
            end
            fooof_data = fooof(freqs, spectra_standing_power(ff,:), f_range, settings, true);
            spectra_standing_fooof(ff,:) = fooof_data.power_spectrum - fooof_data.ap_fit;
            ap_fit_standing(ff,:)        = fooof_data.ap_fit;
        end
        PSD(counter).spectra_standing_fooof_log10 = spectra_standing_fooof;
        PSD(counter).ap_fit_standing_log10        = ap_fit_standing;
        PSD(counter).spectra_standing_fooof_db  = 10*spectra_standing_fooof;
        PSD(counter).ap_fit_standing_db         = 10*ap_fit_standing;
        PSD(counter).ap_fit_standing_power        = 10.^ap_fit_standing;
        PSD(counter).spectra_standing_fooof_power = spectra_standing_power(:,f1:f2) - 10.^ap_fit_standing;
    end
end

filename = fullfile(saveDir,'PSD.mat');
save(filename,'PSD','-v7.3')
disp('done.')