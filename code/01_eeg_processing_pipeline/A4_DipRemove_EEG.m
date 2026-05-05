%% Remove "non-brain" dipoles
%
%  This script performs dipole fitting, if haven't done so already, then
%  flag non-brain dipoles based on PSD upward slope, powpowcat removal 
%  (taken from Jacobsen et al, 2024), outside of brain, and rv < 15, 
%
% The script processes all .set files in the selected 03_bootstrapped_* 
% folder and saves cleaned datasets to 04_dipRemoved_* folder 
% along with a summary table of flagged components.
%
% requires PowPowCAT plugin and Viewprops plugin, and fieldstrip
%
% LM 122225
%%
clear all;
% add libraries
eeglab; ft_defaults; clc; close all;
if ispc
    mainDir     = 'C:\Git\DoD-Gait';
elseif ismac
    mainDir = '/Users/Leo/Git/DoD-Gait';
else
    mainDir = '/home/leo/Documents/DoD-Gait';
end
funcPath    = fullfile(mainDir,'code','func');
addpath(funcPath);

disp('Load in 03_bootstrapped folder.')
figDir      = fullfile(mainDir,'reports','data_quality');
dataDir     = fullfile(mainDir,'data');
workingDir  = uigetdir(dataDir,'Load in 03_bootstrapped folder.');
label       = extractAfter(workingDir,'bootstrapped_');
eegFiles    = dir(fullfile(workingDir,'*.set'));
saveFolder =  fullfile(mainDir,'data',['04_dipRemoved_' label]);
mkdir(saveFolder);

% set up var for tracking in parfor
n = length(eegFiles);
Subject                = strings(n, 1);
Condition              = strings(n, 1);

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

if check % if true, perform mannual check of head model warps
    for s = 1:length(eegFiles)
        if contains(eegFiles(s).name,'rest')
            continue
        end
        EEG = pop_loadset('filename', eegFiles(s).name,  'filepath', eegFiles(s).folder);
        if length(EEG.urchanlocs) == 71
            EEG.urchanlocs(65:end-1) = []; % maintain consistency for STUDY.
        end
        EEG = eeg_checkset(EEG);
        EEG = getHeadModelSettings(EEG); % perform wapring and checking here

        %save EEG
        EEG.setname = [EEG.setname ' - coregistered'];
        pop_saveset(EEG, eegFiles(s).name, eegFiles(s).folder);
    end
else
    for s = 1:length(eegFiles)
        if contains(eegFiles(s).name,'rest')
            continue
        end
        EEG = pop_loadset('filename', eegFiles(s).name,  'filepath', eegFiles(s).folder);     
        if length(EEG.urchanlocs) == 71
            EEG.urchanlocs(65:end-1) = []; % maintain consistency for STUDY.
        end
        EEG = eeg_checkset(EEG);
        
        if isfield(EEG.dipfit,'model')
            disp('Dipfit already performed...skiping.')
            EEG.etc.DigitizedChannels = T(strcmp({T.subject},EEG.subject)).DigitizedChannels;
            EEG.etc.ManuallyFixCoreg = T(strcmp({T.subject},EEG.subject)).ManuallyFixCoreg;
            pop_saveset(EEG, eegFiles(s).name, eegFiles(s).folder);
        else
            disp(['Adding warp transform values for ' EEG.subject]);
            EEG = pop_dipfit_settings(EEG,...
                'coordformat','MNI',...
                'coord_transform',T(strcmp({T.subject},EEG.subject)).warpTransform,...
                'model','standardBEM');
            EEG.etc.DigitizedChannels = T(strcmp({T.subject},EEG.subject)).DigitizedChannels;
            EEG.etc.ManuallyFixCoreg = T(strcmp({T.subject},EEG.subject)).ManuallyFixCoreg;
            EEG = eeg_checkset(EEG);

            %save EEG
            EEG.setname = [EEG.setname ' - coregistered'];
            pop_saveset(EEG, eegFiles(s).name, eegFiles(s).folder);
        end
    end
end

%% run dipfit
delete(gcp('nocreate'))
parpool("Processes",7)
parfor s = 1:length(eegFiles)

    EEG = pop_loadset('filename', eegFiles(s).name,  'filepath', eegFiles(s).folder);
    
    if isempty(EEG.icaact)
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    end

    InitialIC(s) = size(EEG.icaact,1);

    %%% remove dipoles that has an upward slope in PSD
    slopeFlagged = false(1,length(EEG.reject.gcompreject));
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
    varFlagged = false(1,length(EEG.reject.gcompreject)); % not used
    
    %%% perform dipfit
    dipoleXyz       = nan(length(EEG.dipfit.model),3);
    rv_vals         = nan(length(EEG.dipfit.model),1);

    if isfield(EEG.dipfit,'boot')
        disp('Dipfit bootstrap performed... getting bootstrapped xyz and rv values.')
        dipoleXyz   = cell2mat(arrayfun(@(m)m.posxyz(1,:), EEG.dipfit.boot, 'UniformOutput', false)');
        rv_vals     = [EEG.dipfit.boot.rv];
        EEG.etc.bootstrappedUsed = true;
        % note: may return NaN
    elseif isfield(EEG.dipfit,'model') && ~isfield(EEG.dipfit,'boot')
        disp('Dipfit already performed... getting xyz and rv values.')
        dipoleXyz   = cell2mat(arrayfun(@(m)m.posxyz(1,:), EEG.dipfit.model, 'UniformOutput', false)');
        rv_vals     = [EEG.dipfit.model.rv];
        EEG.etc.bootstrappedUsed = false;
    else
        EEG = pop_multifit(EEG,[],'threshold',100,'rmout','on','plotopt',{'normlen','on'});
        dipoleXyz   = cell2mat(arrayfun(@(m)m.posxyz(1,:), EEG.dipfit.model, 'UniformOutput', false)');
        rv_vals     = [EEG.dipfit.model.rv];
        EEG.etc.bootstrappedUsed = false;
    end

    %%% bootstrap flagged %%%
    bootFlagged = [EEG.dipfit.boot.flag];
    
    %%% powpow removal %%%
    powpowFlagged = false(1,length(EEG.reject.gcompreject));
    EEG_powpow = calc_PowPowCAT(EEG,55,2,2,10); % require PowPowCAT plugin
    powpowFlagged_removeIdx = PowPowCat_ICrej(EEG_powpow,1,figDir); % from https://github.com/jacobsen-noelle/ExoAdapt-DualEEG-Processing
    close all;
    powpowFlagged(powpowFlagged_removeIdx) = true;
    
    %%% remove dipoles outside the brain again using ft_sourcedepth %%%
    outsideFlagged = false(length(EEG.reject.gcompreject),1);
    hdm = load(EEG.dipfit.hdmfile); % This returns 'vol'.
    
    % check for nan and fill in outside value for so that ft_ft_sourcedepth
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

    %%% step 2 removal %%%
    totalFlagged = slopeFlagged | varFlagged | powpowFlagged | rvFlagged | outsideFlagged | bootFlagged;
    EEG.reject.gcompreject = totalFlagged;
    EEG.etc.dipRemovedIdx = totalFlagged;
    EEG = pop_subcomp(EEG,'',0,0);          % [] or '' means removing components flagged for rejection
    EEG.dipfit.boot(totalFlagged) = []; % remove ic from bootstrapped data
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
    bootstrap_Flagged(s)        = sum(bootFlagged);
    total_Flagged(s)            = sum(totalFlagged);
    Kept(s)                     = size(EEG.icaact,1); % IC kept
    DigiChans(s)                = EEG.etc.DigitizedChannels;
    ManuallyFixCoreg(s)         = EEG.etc.ManuallyFixCoreg;
    BootstrappedUsed(s)         = EEG.etc.bootstrappedUsed;

    % % save EEG 
    EEG.setname = [EEG.setname ' - DipRemoved'];
    pop_saveset(EEG, EEG.filename, saveFolder);
end
delete(gcp("nocreate"));

%% save table
sT = table(Subject, Condition, DataPercentUsedForICA', DurationUsedForICA', ...
          InitialIC', slope_Flagged', var_Flagged', powpow_Flagged',...
          rv_Flagged', outside_Flagged', bootstrap_Flagged', total_Flagged',...
          Kept', DigiChans', ManuallyFixCoreg', BootstrappedUsed', ...
          'VariableNames', {'Subject', 'Condition', 'DataPercentUsedForICA','DurationUsedForICA',...
          'InitialIC', 'slope_Flagged','var_Flagged','powpow_Flagged', ...
          'rv_Flagged', 'outside_Flagged', 'bootstrap_Flagged', 'total_Flagged',...
          'Kept', 'Digitized_Channels', 'ManuallyFixed', 'BootstrappedUsed'});

writetable(sT,fullfile(figDir,['A4_DipRemoved_' label '.xlsx']));
disp('done.')