%% EEG Dipfit Bootstrap Analysis and Dipole Matching
%
% This script match dipoles across multiple runs and perform bootstrapping
% of dipoles locations and rv. 
% It calculates correlations between ICs from different runs, matches to
% the ref (1st) run.
%
% 1. Load EEG datasets and apply dipole fitting settings.
% 2. Compute correlations between IC activations and topographies across runs.
% 3. Match ICs and identify poorly correlated ICs and mark them as NaN.
% 5. Perform bootstrap resampling to compute robust dipole positions and rv
% 6. Save the processed EEG datasets with bootstrap results.
%
% LM 122325
%%
clear all; clc; close all; warning('off','all');
if ispc
    %mainDir    = 'E:\clab\DoD-Gait';
    mainDir     = 'C:\Git\DoD-Gait';
else
    mainDir    = '/home/leo/Documents/DoD-Gait';
end
eeglab; ft_defaults; clc; close all;

funcPath    = fullfile(mainDir,'code','func');
addpath(funcPath);
figDir      = fullfile(mainDir,'reports','data_quality');
dataDir     = fullfile(mainDir,'data');
AMICAed_folders = dir(fullfile(dataDir,'*AMICAed*'));
% make sure it is in ascending order
nums = cellfun(@(s) str2double(regexp(s,'\d+$','match','once')), {AMICAed_folders.name});
[~, sortidx] = sort(nums,'ascend');
AMICAed_folders = AMICAed_folders(sortidx);

% check folders has all the subjects
for f = 1:length(AMICAed_folders)
    if length(dir(fullfile(AMICAed_folders(1).folder, AMICAed_folders(1).name,'*.set'))) ~= 58
        keyboard
    end
end
% label       = extractAfter(AMICAed_folders(1).name,'AMICAed_');
% label       = extractBefore(label,'_01');
% saveFolder =  fullfile(dataDir,['03_bootstrapped_' label]);
saveFolder =  fullfile(dataDir,'03_bootstrapped');
mkdir(saveFolder)
figDir2 = fullfile(figDir, 'A3_bootstrap_matching');
if ~exist(figDir2, 'dir'), mkdir(figDir2); end

% load in transformation data. 
warpFile = fullfile(dataDir,'Subject_warpTransformData.mat');
load(warpFile); % struct is called T.

% parameters
corrThres = 0.9; % matching dipoles based on icaact and icawinv

%% co-register
check = true;
eegFiles = dir(fullfile(AMICAed_folders(1).folder, AMICAed_folders(1).name,'*.set'));
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
    T = struct();
    for s = 1:length(eegFiles)
        EEG = pop_loadset('filename', eegFiles(s).name,  'filepath', eegFiles(s).folder);
        if length(EEG.urchanlocs) == 71
            EEG.urchanlocs(65:end-1) = []; % maintain consistency for STUDY.
        end
        EEG = eeg_checkset(EEG);
        EEG = getHeadModelSettings(EEG); % perform wapring and checking here

        % save EEG
        %EEG.setname = [EEG.setname ' - coregistered'];
        %pop_saveset(EEG, eegFiles(s).name, eegFiles(s).folder);

        % save coordinate transform values
        T(s).subject                  = EEG.subject;
        T(s).group                    = EEG.group;
        T(s).warpTransform            = EEG.dipfit.coord_transform;
    end
    %%% save T 'Subject_warpTransformData.mat' here
    keyboard
else
    for s = 1:length(eegFiles)
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

%% main loop
delete(gcp('nocreate'))
parpool("Processes",8)
parfor sub = 1:length(T)
    boot = [];
    subject = T(sub).subject;

    dip = struct('run', [], ...
        'ref_idx', [], ...
        'curr_idx', [], ...
        'posxyz', [], ...
        'momxyz', [], ...
        'rv', [], ...
        'act_corr', [], ...
        'topo_corr', [],...
        'icawinv', [],...
        'poor_correlation_idx', []);

    for run = 1:length(AMICAed_folders)
        label2       = extractAfter(AMICAed_folders(run).name,'AMICAed_');
        disp(['...running ' label2]);
        workingDir = fullfile(AMICAed_folders(run).folder,AMICAed_folders(run).name);
        eegFile    = dir(fullfile(workingDir,['*' subject '*.set']));

        EEG = pop_loadset('filename', eegFile.name,  'filepath', eegFile.folder);
        disp(['Adding warp transform values for ' subject]);
        EEG = pop_dipfit_settings(EEG,...
            'coordformat','MNI',...
            'coord_transform',T(strcmp({T.subject},EEG.subject)).warpTransform,...
            'model','standardBEM');
        EEG.etc.DigitizedChannels = T(strcmp({T.subject},EEG.subject)).DigitizedChannels;
        EEG.etc.ManuallyFixCoreg = T(strcmp({T.subject},EEG.subject)).ManuallyFixCoreg;
        EEG = eeg_checkset(EEG);
        EEG = pop_multifit(EEG,[],'threshold',100,'rmout','on','plotopt',{'normlen','on'});

        if run == 1
            ref_dipfit  = EEG.dipfit; % save dipfitted from ref set
            ref_icaact  = EEG.icaact;
            ref_icawinv = EEG.icawinv;
            for idx = 1:length(EEG.dipfit.model)
                dip(run).run = run;
                dip(run).ref_idx(idx) = idx;
                dip(run).curr_idx(idx) = idx;
                dip(run).posxyz(idx,:) = EEG.dipfit.model(idx).posxyz;
                dip(run).momxyz(idx,:) = EEG.dipfit.model(idx).momxyz;
                dip(run).rv(idx) = EEG.dipfit.model(idx).rv;
                dip(run).act_corr(idx) = nan;
                dip(run).topo_corr(idx) = nan;
                dip(run).icawinv(:,idx) = ref_icawinv(:,idx);
            end
        else
            curr_icaact   = EEG.icaact;
            curr_icawinv  = EEG.icawinv;

            % get correlation array
            ref_z      = zscore(ref_icaact,0,2); % ignore amp but maintain shape
            curr_z    = zscore(curr_icaact,0,2);
            corrMatrix  = corr(ref_z', curr_z'); % get correlation of all matches
            num_ref    = size(ref_z, 1);
            num_curr  = size(curr_z, 1);

            %% Matchpairs method
            costMatrix = max(abs(corrMatrix(:))) - abs(corrMatrix);
            costUnmatched = max(costMatrix(:)) * 10;
            [M, uR, uC]   = matchpairs(costMatrix, costUnmatched, 'min'); % row = ref, col = curr
            matchedCorr   = arrayfun(@(i) corrMatrix(M(i,1), M(i,2)), 1:size(M,1))';
            matchedCorr     = abs(matchedCorr);

            currIdxMap = nan(1,num_ref);
            currIdxMap(M(:,1)) = M(:,2);

            for refIdx = 1:num_ref
                currIdx = currIdxMap(refIdx);

                fprintf('Ref %d ↔ Current %d, corr=%.3f\n', ...
                    refIdx, currIdx, corrMatrix(refIdx, currIdx));
                dip(run).run                    = run;
                dip(run).ref_idx(refIdx)        = refIdx;
                dip(run).curr_idx(refIdx)       = currIdx;
                dip(run).act_corr(refIdx)       = abs(corrMatrix(refIdx, currIdx));
                dip(run).topo_corr(refIdx)      = abs(corr(ref_icawinv(:,refIdx), ...
                    curr_icawinv(:,currIdx)));
                dip(run).posxyz(refIdx,:)       = EEG.dipfit.model(currIdx).posxyz;
                dip(run).momxyz(refIdx,:)       = EEG.dipfit.model(currIdx).momxyz;
                dip(run).rv(refIdx)             = EEG.dipfit.model(currIdx).rv;
                dip(run).icawinv(:,refIdx)      = curr_icawinv(:,currIdx);
            end

            %%% plot correlation
            fig = figure('Color','w','WindowState', 'maximized','Visible','off');
            nexttile
            imagesc(corrMatrix);
            colorbar;
            title('ICA act correlation between runs');
            xlabel('Current run ICs');
            ylabel('REF run ICs');

            nexttile
            x = dip(run).act_corr;
            y = dip(run).topo_corr;

            % Indices for good vs low correlations
            goodIdx = x >= corrThres | y >= corrThres;
            badIdx  = x < corrThres &  y < corrThres;

            dip(run).poor_correlation_idx = badIdx;

            scatter(x(goodIdx), y(goodIdx), 'b', 'filled', 'SizeData',50);hold on
            scatter(x(badIdx), y(badIdx), 'r', 'filled', 'SizeData',75);
            xlabel('Activation correlation');
            ylabel('Topography correlation');
            title('Correlation per IC');
            grid on;
            numICs = length(dip(run).act_corr);
            for i = 1:numICs
                text(dip(run).act_corr(i), dip(run).topo_corr(i), ...
                    sprintf('%d', i), ...
                    'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'right', ...
                    'FontSize', 8, 'Color', 'k');
            end

            if sum(badIdx) > 0
                figName = fullfile(figDir2, sprintf('FLAGGED_%s_corr_run-%d.png', EEG.subject, run));
            else
                figName = fullfile(figDir2, sprintf('%s_corr_run-%d.png', EEG.subject, run));
            end

            exportgraphics(fig, figName, 'Resolution', 600);
            sgtitle(sprintf('%s_corr_run-%d', EEG.subject, run), 'Interpreter', 'none')
            close(fig)
        end
    end

    ic_num = length(dip(1).rv);
    pos_x_vals = []; pos_y_vals = []; pos_z_vals = [];
    mom_x_vals = []; mom_y_vals = []; mom_z_vals = [];

    for rr = 1:length(dip)
        pos_x_vals(rr,:) = dip(rr).posxyz(:,1);
        pos_y_vals(rr,:) = dip(rr).posxyz(:,2);
        pos_z_vals(rr,:) = dip(rr).posxyz(:,3);
        mom_x_vals(rr,:) = dip(rr).momxyz(:,1);
        mom_y_vals(rr,:) = dip(rr).momxyz(:,2);
        mom_z_vals(rr,:) = dip(rr).momxyz(:,3);
    end
    rv_vals = cat(1, dip.rv);
    icawinv_vals = cat(3, dip.icawinv);

    poor_corr = cat(1, dip.poor_correlation_idx);
    poor_corr = [zeros(1,size(poor_corr,2)); poor_corr]; % add first/ref row
    poor_corr_sum = sum(poor_corr);
    poor_corr = logical(poor_corr); % convert to logical

    % mark nan for poorly correlated components
    pos_x_vals(poor_corr)   = nan;
    pos_y_vals(poor_corr)   = nan;
    pos_z_vals(poor_corr)   = nan;
    mom_x_vals(poor_corr)   = nan;
    mom_y_vals(poor_corr)   = nan;
    mom_z_vals(poor_corr)   = nan;

    rv_vals(poor_corr)  = nan;
    nChan = size(icawinv_vals,1);
    nanMask = permute(isnan(poor_corr),[3 2 1]);
    nanMask = repmat(nanMask,[nChan 1 1]);
    icawinv_vals(nanMask) = NaN; % [chan #, ic #, run #]

    for ic = 1:ic_num

        pos_x_temp  = pos_x_vals(:,ic);
        pos_y_temp  = pos_y_vals(:,ic);
        pos_z_temp  = pos_z_vals(:,ic);

        mom_x_temp  = mom_x_vals(:,ic);
        mom_y_temp  = mom_y_vals(:,ic);
        mom_z_temp  = mom_z_vals(:,ic);

        rv_temp = rv_vals(:,ic);
        icawinv_temp = squeeze(icawinv_vals(:,ic,:));

        nboot = 2000;
        N = length(dip);
        boot_pos_x  = nan(1,nboot);
        boot_pos_y  = nan(1,nboot);
        boot_pos_z  = nan(1,nboot);
        boot_mom_x  = nan(1,nboot);
        boot_mom_y  = nan(1,nboot);
        boot_mom_z  = nan(1,nboot);
        boot_rv     = nan(1,nboot);
        boot_icawinv = nan(size(icawinv_temp,1), nboot);

        % make sure icawinv is not flipped
        first_map = icawinv_temp(:,1);
        for r = 1:size(icawinv_temp,2)   % loop over runs
            curr_map = icawinv_temp(:,r);

            % Compute correlation with reference
            if corr(curr_map, first_map) < 0
                icawinv_temp(:,r) = -curr_map;  % flip sign
            end
        end

        for b = 1:nboot
            rng(5489, 'mt19937ar');
            bootIdx = randsample(N, N, true); % resample with replacement
            boot_pos_x(b)  = mean(pos_x_temp(bootIdx),'omitmissing');
            boot_pos_y(b)  = mean(pos_y_temp(bootIdx),'omitmissing');
            boot_pos_z(b)  = mean(pos_z_temp(bootIdx),'omitmissing');
            boot_mom_x(b)  = mean(mom_x_temp(bootIdx),'omitmissing');
            boot_mom_y(b)  = mean(mom_y_temp(bootIdx),'omitmissing');
            boot_mom_z(b)  = mean(mom_z_temp(bootIdx),'omitmissing');
            boot_rv(b)      = mean(rv_temp(bootIdx),'omitmissing');
            boot_icawinv(:,b) = mean(icawinv_temp(:,bootIdx),2,'omitmissing');
        end

        boot(ic).posxyz  = [median(boot_pos_x,'omitmissing') median(boot_pos_y,'omitmissing') median(boot_pos_z,'omitmissing')];
        boot(ic).momxyz  = [median(boot_mom_x,'omitmissing') median(boot_mom_y,'omitmissing') median(boot_mom_z,'omitmissing')];
        boot(ic).rv = median(boot_rv,'omitmissing');
        boot(ic).icawinv = median(boot_icawinv,2,'omitmissing');

        boot(ic).ci_pos_x       = prctile(boot_pos_x, [2.5 97.5]);
        boot(ic).ci_pos_y       = prctile(boot_pos_y, [2.5 97.5]);
        boot(ic).ci_pos_z       = prctile(boot_pos_z, [2.5 97.5]);
        boot(ic).ci_mom_x       = prctile(boot_mom_x, [2.5 97.5]);
        boot(ic).ci_mom_y       = prctile(boot_mom_y, [2.5 97.5]);
        boot(ic).ci_mom_z       = prctile(boot_mom_z, [2.5 97.5]);
        boot(ic).ci_rv      = prctile(boot_rv, [2.5 97.5]);
        boot(ic).ci_icawinv = prctile(boot_icawinv, [2.5 97.5], 2);

        boot(ic).flag = false;
        if sum(isnan(rv_temp)) > 1
            boot(ic).flag = true; % flag ic with more than 1 poorly correlated runs
        end
    end

    % add bootstrap values to EEG.dipfit
    folder1 = [extractBefore(eegFile.folder,'moreAffected') 'moreAffected_01'];
    EEG = pop_loadset('filename', eegFile.name,  'filepath', folder1);
    EEG.dipfit = ref_dipfit;
    EEG.dipfit.boot = boot;
    EEG.etc.bootstrap_dipfit_info = dip;
    setname = extractBefore(EEG.filename,'.set');
    EEG.setname = [setname ' - dipfitted - bootstrapped'];
    pop_saveset(EEG, EEG.filename, saveFolder);

    %% plot correlated icaact %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figDir3 = fullfile(figDir, 'A3_bootstrap_matching', EEG.subject);
    if ~exist(figDir3, 'dir'), mkdir(figDir3); end

    for rr = 2:length(dip)

        for comp = 1:length(dip(rr).curr_idx)
            disp(['...plotting matching comp ' num2str(comp)]);
            fig = figure('Name', sprintf('Comp %d', comp), ...
                'NumberTitle', 'off', ...
                'WindowState', 'maximized',...
                'Color','w','Visible', 'on');

            % ref topoplot
            ax = nexttile;
            topoplot(EEG.icawinv(:,comp), EEG.chanlocs, 'electrodes', 'off');
            axis on
            keyboard
            xlabel(ax,sprintf(['XYZ: %s \n' ...
                'rv: %.2f'], ...
                sprintf('%.2f ', EEG.dipfit.model(comp).posxyz), ...
                EEG.dipfit.model(comp).rv), ...
                'FontSize', 14);
            title(sprintf('Ref Comp: %2.0f', ...
                comp),'FontSize',22);

            % current topoplot
            ax2 = nexttile;
            topoplot(dip(rr).icawinv(:,comp), EEG.chanlocs, 'electrodes', 'off');
            axis on
            xlabel(ax2,sprintf(['XYZ: %s \n' ...
                'rv: %.2f'], ...
                sprintf('%.2f ', [dip(rr).posxyz(comp,:)]),...
                dip(rr).rv(comp)), ...
                'FontSize', 14);
            title(sprintf('Curr run comp: %2.0f', ...
                dip(rr).curr_idx(comp)),'FontSize',22);

            sgtitle(['Run - ' num2str(rr)],'FontSize',24);

            set(gcf,'Color','w','WindowState', 'maximized');
            % save plot
            if ~boot(comp).flag
                figName = fullfile(figDir3, sprintf('%s_Bootstrap-matching-GOOD_run-%d_comp-%d.png', EEG.subject, rr, comp));
                
                % saveas(fig, figName); % commented to save space
                if strcmp(EEG.subject,'sub-PD04')
                    saveas(fig, figName);
                end
                
                close(fig);
            else
                figName = fullfile(figDir3, sprintf('%s_Bootstrap-matching-BAD_run-%d_comp-%d.png', EEG.subject, rr, comp));
                saveas(fig, figName);
                close(fig);
            end

        end
    end
end
delete(gcp('nocreate'))
disp('done.')