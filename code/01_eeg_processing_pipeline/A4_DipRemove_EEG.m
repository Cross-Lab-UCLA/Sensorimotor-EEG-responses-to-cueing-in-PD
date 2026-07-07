%% Remove "non-brain" dipoles
%
%  This script performs dipole fitting, if haven't done so already, then
%  flag non-brain dipoles based on PSD upward slope, powpowcat removal 
%  (taken from Jacobsen et al, 2024), outside of brain, and rv > 15%.
%
% The script processes all .set files in the selected 03_cluster_* 
% folder and saves cleaned datasets to 04_dipRemoved_* folder 
% along with a summary table of flagged components.
%
% requires PowPowCAT plugin and Viewprops plugin, and fieldstrip
%
% LM 071426
% - moved away from the bootstrap approach from A3
%
%%
clear all;
% add libraries
eeglab; ft_defaults; clc; close all;
if ispc
    mainDir = 'C:\clab\DoD-Gait';
elseif ismac
    mainDir = '/Users/Leo/Git/DoD-Gait';
else
    mainDir = '/home/leo/Documents/DoD-Gait';
end
funcPath    = fullfile(mainDir,'code','01_eeg_processing_pipeline','func');
addpath(funcPath);

disp('Load in 03_cluster folder.')
figDir      = fullfile(mainDir,'reports','data_quality');
dataDir     = fullfile(mainDir,'data');
workingDir  = uigetdir(dataDir,'Load in 03_cluster folder.');
eegFiles    = dir(fullfile(workingDir,'*.set'));

[~, folderName] = fileparts(workingDir);
label       = extractAfter(workingDir,'cluster_');
saveFolder = fullfile(mainDir, 'data', ['04_dipRemoved', repmat(['_' label], 1, ~isempty(label))]);
mkdir(saveFolder);

% set up var for tracking in parfor
n = length(eegFiles);
Subject                = strings(n, 1);
Condition              = strings(n, 1);
InitialIC              = NaN(n, 1);
DataPercentUsedForICA  = NaN(n, 1);
DurationUsedForICA     = NaN(n, 1);
slope_Flagged          = NaN(n, 1);
var_Flagged            = NaN(n, 1);
powpow_Flagged         = NaN(n, 1);
rv_Flagged             = NaN(n, 1);
outside_Flagged        = NaN(n, 1);
medianOutside_Flagged  = NaN(n, 1);
medianRv_Flagged       = NaN(n, 1);
minValidRuns           = NaN(n, 1);   % fewest valid-dipole runs among a subject's kept clusters
cluster_Flagged        = NaN(n, 1);
total_Flagged          = NaN(n, 1);
Kept                   = NaN(n, 1);
DigiChans              = false(n, 1);
ManuallyFixCoreg       = false(n, 1);
BootstrappedUsed       = false(n, 1);

%% co-register
check = true;
EEG = pop_loadset('filename', eegFiles(end).name,  'filepath', eegFiles(end).folder);
warpFile = fullfile(dataDir,'Subject_warpTransformData.mat');

if isfile(warpFile)
    choice = questdlg('Co-registration warp values are found. Would you like to use them?', ...
        'Select', ...
        'Yes','No','No');
    switch choice
        case 'Yes'
            disp('Leading previously used co-registrations');
            load(warpFile); % struct is called T.
            check = false;
        case 'No'
            disp('Performing co-registration check');
    end
end

if check % if true, perform manual check of head model warps
    T = struct();
    for s = 1:length(eegFiles)
        EEG = pop_loadset('filename', eegFiles(s).name,  'filepath', eegFiles(s).folder);
        if length(EEG.urchanlocs) == 71
            EEG.urchanlocs(65:end-1) = []; % maintain consistency for STUDY.
        end
        EEG = eeg_checkset(EEG);
        disp(['Performing coregister on ' EEG.subject]);
        EEG = getHeadModelSettings(EEG); % perform warping and checking here

        % save EEG
        %EEG.setname = [EEG.setname ' - coregistered'];
        %pop_saveset(EEG, eegFiles(s).name, eegFiles(s).folder);

        % save coordinate transform values
        T(s).subject                  = EEG.subject;
        T(s).group                    = EEG.group;
        T(s).warpTransform            = EEG.dipfit.coord_transform;
        T(s).DigitizedChannels        = EEG.etc.DigitizedChannels;
        T(s).ManuallyFixCoreg            = EEG.etc.ManuallyFixCoreg;

    end
    %%% save T 'Subject_warpTransformData.mat' here
    disp('Saving Warp Transform Table.')
    save(fullfile(dataDir, 'Subject_warpTransformData.mat'), 'T');

else
    % add coord_transform from table
    for s = 1:length(eegFiles)
        EEG = pop_loadset('filename', eegFiles(s).name,  'filepath', eegFiles(s).folder);
        if length(EEG.urchanlocs) == 71
            EEG.urchanlocs(65:end-1) = []; % maintain consistency for STUDY.
        end
        EEG = eeg_checkset(EEG);

        subIdx = strcmp({T.subject}, EEG.subject);
        if ~any(subIdx)
            error('Subject %s not found in table T', EEG.subject);
        end

        if isfield(EEG,'dipfit') && isfield(EEG.dipfit,'coord_transform') ...
                && isequal(EEG.dipfit.coord_transform, T(subIdx).warpTransform)
            disp('coord_transform is already imported, skipping...')
            continue
        end

        % set up dipfit
        EEG = pop_dipfit_settings(EEG,...
            'coordformat','MNI',...
            'coord_transform',T(subIdx).warpTransform,...
            'model','standardBEM');

        EEG.etc.DigitizedChannels = T(subIdx).DigitizedChannels;
        EEG.etc.ManuallyFixCoreg = T(subIdx).ManuallyFixCoreg;

        %save EEG
        EEG = eeg_checkset(EEG);
        EEG.setname = [EEG.setname ' - coregistered'];
        pop_saveset(EEG, eegFiles(s).name, eegFiles(s).folder);
    end

end

%% run removal processing
delete(gcp('nocreate'))
parpool("Processes",7)
parfor s = 1:length(eegFiles)

    EEG = pop_loadset('filename', eegFiles(s).name,  'filepath', eegFiles(s).folder);
    
    if isempty(EEG.icaact)
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    end

    InitialIC(s) = size(EEG.icaact,1);
    nIC = size(EEG.icaact,1);   % component count; all flag arrays sized to this

    %%% remove dipoles that has an upward slope in PSD
    slopeFlagged = false(1,nIC);
    for di = 1:size(EEG.icaact,1)
        inputData = EEG.icaact(di,:);
        [spectra1,freqs] = spectopo(inputData, 0, EEG.srate,...
            'winsize', EEG.srate, 'nfft', EEG.srate,...
            'overlap', .5*EEG.srate, 'plot', 'off');
        freqs1 = find(freqs == 4);
        freqs2 = find(freqs == 40);
        spectraCrop = spectra1(freqs1:freqs2);
        p = polyfit(freqs1:freqs2, spectraCrop, 1);
        if p(1) >= 0 % if the dipole has a positive slope, it is removed
            slopeFlagged(di) = true;
        end

        % nexttile
        % plot(freqs1:freqs2,spectraCrop); hold on
        % y3 = polyval(p,freqs1:freqs2);
        % plot(freqs1:freqs2,y3)
    end

    %%% remove ic that is explained more than 90% by a single channel 
    % mix = cat(2,EEG.dipfit.boot.icawinv);
    % chanvar_components = max(mix.^2 ./ sum(mix.^2,1),[],1)';
    % bad_var_components = chanvar_components > 0.9;
    % varFlagged = false(1,size(mix,2));
    % varFlagged(bad_var_components) = true;
    varFlagged = false(1,nIC); % not used
    
    %%% perform dipfit
    if isfield(EEG.dipfit,'cltr')
        disp('Dipfit clustering performed... getting representive xyz and rv values.')
        dipoleXyz   = cell2mat(arrayfun(@(m)m.posxyz(1,:), EEG.dipfit.cltr, 'UniformOutput', false)');
        rv_vals     = [EEG.dipfit.cltr.rv];

        nClu = length(EEG.dipfit.cltr);
        dipoleXyz_runs = nan(nClu,3);
        rv_runs        = nan(nClu,1);
        validRuns      = nan(nClu,1);
        for k = 1:nClu
            pr = EEG.dipfit.cltr(k).pos_runs;
            dipoleXyz_runs(k,:) = median(pr,'omitnan'); % median over valid runs
            rv_runs(k)          = median(EEG.dipfit.cltr(k).rv_runs,'omitnan');
            validRuns(k)        = sum(all(~isnan(pr),2)); % check for NaNs
            EEG.dipfit.cltr(k).nValidRuns = validRuns(k);
        end
        EEG.etc.clusterUsed = true;

    elseif isfield(EEG.dipfit,'model') && ~isfield(EEG.dipfit,'cltr')
        disp('Dipfit already performed... getting xyz and rv values.')
        dipoleXyz   = cell2mat(arrayfun(@(m)m.posxyz(1,:), EEG.dipfit.model, 'UniformOutput', false)');
        rv_vals     = [EEG.dipfit.model.rv];
        EEG.etc.clusterUsed = false;
    else
        EEG = pop_multifit(EEG,[],'threshold',100,'rmout','on','plotopt',{'normlen','on'});
        dipoleXyz   = cell2mat(arrayfun(@(m)m.posxyz(1,:), EEG.dipfit.model, 'UniformOutput', false)');
        rv_vals     = [EEG.dipfit.model.rv];
        EEG.etc.clusterUsed = false;
    end

    %%% bootstrap flagged %%%
    % if isfield(EEG.dipfit,'cltr')
    %     bootFlagged = [EEG.dipfit.cltr.flag];
    % else
    %     bootFlagged = false(1,length(EEG.reject.gcompreject));
    % end
    
    %%% powpow removal %%%
    powpowFlagged = false(1,nIC);
    EEG_powpow = calc_PowPowCAT(EEG,55,2,2,10); % require PowPowCAT plugin
    powpowFlagged_removeIdx = PowPowCat_ICrej(EEG_powpow,1,figDir); % from https://github.com/jacobsen-noelle/ExoAdapt-DualEEG-Processing
    close all;
    powpowFlagged(powpowFlagged_removeIdx) = true;
    
    %%% remove dipoles outside the brain again using ft_sourcedepth %%%
    outsideFlagged = false(nIC,1);
    hdm = load(EEG.dipfit.hdmfile); % This returns 'vol'.
    
    % check for nan and fill in outside value for so that ft_sourcedepth
    % doesn't crash
    nanIdx = find(isnan(dipoleXyz(:,1)));
    if ~isempty(nanIdx)
        dipoleXyz(nanIdx,:) = repmat([99 99 99], numel(nanIdx), 1);
    end

    depth = round(ft_sourcedepth(dipoleXyz, hdm.vol),1);
    % from ft_sourcedepth
    %   A negative depth indicates that the source is inside the source
    %   compartment, positive indicates outside.
    depthThreshold = 0;
    insideBrainIdx = find(depth < depthThreshold);
    outsideFlagged = ~ismember(1:size(EEG.icaact,1), insideBrainIdx);

    %%% remove RV > .15 %%%
    rvThreshold = .15;
    rvFlagged = rv_vals > rvThreshold | isnan(rv_vals);

    %%% median-across-runs requirement (clusters only) %%%
    % A cluster passes only if its per-run MEDIAN dipole is also inside the
    % brain and has median RV < threshold.
    medianOutsideFlagged = false(1,nIC);
    medianRvFlagged      = false(1,nIC);
    if EEG.etc.clusterUsed
        medXyz = dipoleXyz_runs;                       % [nIC x 3] per-run median position
        nanIdxM = find(isnan(medXyz(:,1)));            % all-NaN clusters -> force outside
        if ~isempty(nanIdxM)
            medXyz(nanIdxM,:) = repmat([99 99 99], numel(nanIdxM), 1);
        end
        depthM   = round(ft_sourcedepth(medXyz, hdm.vol),1);
        insideM  = find(depthM < depthThreshold);
        medianOutsideFlagged = ~ismember(1:nIC, insideM);          % outside brain (row)
        medianRvFlagged      = rv_runs' > rvThreshold | isnan(rv_runs');
    end

    %%% step 2 removal %%%
    % totalFlagged = slopeFlagged | varFlagged | powpowFlagged | rvFlagged | outsideFlagged ...
    %     | medianOutsideFlagged | medianRvFlagged;    % including cluster medians
    
    totalFlagged = slopeFlagged | varFlagged | powpowFlagged | rvFlagged | outsideFlagged;

    EEG.reject.gcompreject = totalFlagged;
    EEG.etc.dipRemovedIdx = totalFlagged;
    EEG = pop_subcomp(EEG,'',0,0);          % [] or '' means removing components flagged for rejection
    
    EEG.icaweights = pinv(EEG.icawinv);
    % A3 builds each IC from a different AMICA run, so pop_subcomp can leak part of each
    % kept IC into the others. Leakage is negligible. Nonetheless, unmixing with pinv of 
    % the kept scalp maps undoes this.
    % Alternative: skip pop_subcomp and prune icaweights/icawinv/dipfit.model/cltr with
    % ~totalFlagged directly; channel data then stay uncleaned, but the ICs are unchanged.
    
    if isfield(EEG.dipfit,'cltr')
        EEG.dipfit.cltr(totalFlagged) = []; % remove ic from IC cluster data
    end
    EEG = eeg_checkset(EEG);

    % write to table
    Subject(s)                  = EEG.subject;
    Condition(s)                = EEG.condition;
    DataPercentUsedForICA(s)    = EEG.etc.DataPercentUsedForICA;
    DurationUsedForICA(s)       = EEG.etc.DurationUsedForICA;
    slope_Flagged(s)            = sum(slopeFlagged);
    var_Flagged(s)              = sum(varFlagged);
    powpow_Flagged(s)           = sum(powpowFlagged);
    rv_Flagged(s)               = sum(rvFlagged);
    outside_Flagged(s)          = sum(outsideFlagged);
    medianOutside_Flagged(s)    = sum(medianOutsideFlagged);
    medianRv_Flagged(s)         = sum(medianRvFlagged);
    if EEG.etc.clusterUsed && any(~totalFlagged)
        minValidRuns(s)         = min(validRuns(~totalFlagged));   % weakest kept cluster
    end
    %cluster_Flagged(s)          = sum(clusterFlagged);
    total_Flagged(s)            = sum(totalFlagged);
    Kept(s)                     = size(EEG.icaact,1); % IC kept
    DigiChans(s)                = EEG.etc.DigitizedChannels;
    ManuallyFixCoreg(s)         = EEG.etc.ManuallyFixCoreg;
    BootstrappedUsed(s)         = EEG.etc.clusterUsed;

    % % save EEG 
    EEG.setname = [EEG.setname ' - DipRemoved'];
    pop_saveset(EEG, EEG.filename, saveFolder);
end
delete(gcp("nocreate"));

%% save table
sT = table(Subject, Condition, DataPercentUsedForICA, DurationUsedForICA, ...
          InitialIC, slope_Flagged, var_Flagged, powpow_Flagged,...
          outside_Flagged, medianOutside_Flagged, rv_Flagged, medianRv_Flagged, ...
          minValidRuns, cluster_Flagged, total_Flagged,...
          Kept, DigiChans, ManuallyFixCoreg, BootstrappedUsed, ...
          'VariableNames', {'Subject', 'Condition', 'DataPercentUsedForICA','DurationUsedForICA',...
          'InitialIC', 'slope_Flagged','var_Flagged','powpow_Flagged', ...
          'outside_Flagged', 'medianOutside_Flagged(notUsed)', 'rv_Flagged', 'medianRv_Flagged(notUsed)', ...
          'minValidRuns', 'cluster_Flagged', 'total_Flagged',...
          'Kept', 'Digitized_Channels', 'ManuallyFixed', 'BootstrappedUsed'});

writetable(sT,fullfile(figDir,'A4_DipRemoved.xlsx'));
disp('done.')