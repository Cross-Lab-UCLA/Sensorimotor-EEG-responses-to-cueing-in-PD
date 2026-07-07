%% Make cued averaged group PSD - cluster-based permutation throughout
% This script performs PSD analysis and plot.
%
% Updated A12_PSD_Stats_and_Plot.m that runs FieldTrip cluster-based
% permutation statistics (over the frequency dimension) for EVERY effect,
% and tests the group x cue interaction with per-subject difference scores.
% See: https://www.fieldtriptoolbox.org/faq/stats/clusterstats_interaction/
%
% Testing using statcondfieldtrip permutation + 'mcorrect','cluster':
%   Group main   : indepsamplesF over 3 groups (cue-averaged)
%   Group pairwise (x3)      : indepsamplesT (cue-averaged)
%   Cue main     : depsamplesFunivariate over 3 cues (group-pooled)
%   Cue pairwise (x3)        : depsamplesT (group-pooled)
%   Group x cue interaction  : per-subject cue-contrast difference spectra,
%                              then between-groups tests on those diffs:
%       audi-nocue, visu-nocue, visu-audi  x  {omnibus indepsamplesF,
%                                              pairwise indepsamplesT}
%
% LM 093026
%%
clear all; clc; eeglab; close all
rng(134);
if ispc
    mainDir     = 'C:\clab\DoD-Gait';
elseif ismac
    mainDir     = '/Users/Leo/Git/DoD-Gait';
else
    mainDir     = '/home/leo/Documents/DoD-Gait';
end
funcPath    = fullfile(mainDir,'code','01_eeg_processing_pipeline','func');
addpath(funcPath)

% directories
dataDir = fullfile(mainDir, 'data');
figFolder = fullfile(mainDir,'results','figures');

% load demo data
T_demo = readtable(fullfile(dataDir,"demo_data_analysis.xlsx"));

% load data
psdFile  = dir(fullfile(dataDir, '08_PSD','*.mat'));
load(fullfile(psdFile.folder,psdFile.name));

% stats params
naccu   = 20000;
alpha   = 0.05;

% false plots every significant post-hoc, including under a null omnibus
gatePostHocOnOmnibus = true;

%% organize
clusterLabel = {'LSM' 'RSM'};
PSD_organized.LSM = [];
PSD_organized.RSM = [];

for clu = 1:length(clusterLabel)
    idx = find(strcmp({PSD.roi},clusterLabel{clu}));
    PSDtmp = PSD(idx);
    hc_idx = find(strcmp({PSDtmp.subgroup},'hc'));
    fog_idx = find(strcmp({PSDtmp.subgroup},'fog'));
    nofog_idx = find(strcmp({PSDtmp.subgroup},'nofog'));

    %%% HC %%%
    PSD_db_hc = PSDtmp(hc_idx);
    PSD_db_hc_nocue_walking     = vertcat(arrayfun(@(s) s.spectra_fooof_db(1, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_nocue_walking     = vertcat(PSD_db_hc_nocue_walking{:});
    PSD_db_hc_nocue_standing    = vertcat(arrayfun(@(s) s.spectra_standing_fooof_db(1, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_nocue_standing    = vertcat(PSD_db_hc_nocue_standing{:});
    PSD_db_hc_nocue             = PSD_db_hc_nocue_walking - PSD_db_hc_nocue_standing;

    PSD_db_hc_audi_walking  = vertcat(arrayfun(@(s) s.spectra_fooof_db(2, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_audi_walking  = vertcat(PSD_db_hc_audi_walking{:});
    PSD_db_hc_audi_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof_db(2, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_audi_standing = vertcat(PSD_db_hc_audi_standing{:});
    PSD_db_hc_audi          = PSD_db_hc_audi_walking - PSD_db_hc_audi_standing;

    PSD_db_hc_visu_walking  = vertcat(arrayfun(@(s) s.spectra_fooof_db(3, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_visu_walking  = vertcat(PSD_db_hc_visu_walking{:});
    PSD_db_hc_visu_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof_db(3, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_visu_standing = vertcat(PSD_db_hc_visu_standing{:});
    PSD_db_hc_visu          = PSD_db_hc_visu_walking-PSD_db_hc_visu_standing;

    PSD_organized.(clusterLabel{clu}).subjects_hc = {PSDtmp(hc_idx).subject}';
    PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,1) = PSD_db_hc_nocue;
    PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,2) = PSD_db_hc_audi;
    PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,3) = PSD_db_hc_visu;

    PSD_db_hc_nocue_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit_db(1, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_nocue_ap_fit  = vertcat(PSD_db_hc_nocue_ap_fit{:});
    PSD_db_hc_audi_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit_db(2, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_audi_ap_fit  = vertcat(PSD_db_hc_audi_ap_fit{:});
    PSD_db_hc_visu_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit_db(3, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_visu_ap_fit  = vertcat(PSD_db_hc_visu_ap_fit{:});

    PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit(:,:,1) = PSD_db_hc_nocue_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit(:,:,2) = PSD_db_hc_audi_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit(:,:,3) = PSD_db_hc_visu_ap_fit;

    %%% nofog %%%
    PSD_db_nofog = PSDtmp(nofog_idx);
    PSD_db_nofog_nocue_walking     = vertcat(arrayfun(@(s) s.spectra_fooof_db(1, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_nocue_walking     = vertcat(PSD_db_nofog_nocue_walking{:});
    PSD_db_nofog_nocue_standing    = vertcat(arrayfun(@(s) s.spectra_standing_fooof_db(1, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_nocue_standing    = vertcat(PSD_db_nofog_nocue_standing{:});
    PSD_db_nofog_nocue             = PSD_db_nofog_nocue_walking - PSD_db_nofog_nocue_standing;

    PSD_db_nofog_audi_walking  = vertcat(arrayfun(@(s) s.spectra_fooof_db(2, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_audi_walking  = vertcat(PSD_db_nofog_audi_walking{:});
    PSD_db_nofog_audi_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof_db(2, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_audi_standing = vertcat(PSD_db_nofog_audi_standing{:});
    PSD_db_nofog_audi          = PSD_db_nofog_audi_walking - PSD_db_nofog_audi_standing;

    PSD_db_nofog_visu_walking  = vertcat(arrayfun(@(s) s.spectra_fooof_db(3, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_visu_walking  = vertcat(PSD_db_nofog_visu_walking{:});
    PSD_db_nofog_visu_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof_db(3, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_visu_standing = vertcat(PSD_db_nofog_visu_standing{:});
    PSD_db_nofog_visu          = PSD_db_nofog_visu_walking - PSD_db_nofog_visu_standing;

    PSD_organized.(clusterLabel{clu}).subjects_nofog = {PSDtmp(nofog_idx).subject}';
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,1) = PSD_db_nofog_nocue;
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,2) = PSD_db_nofog_audi;
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,3) = PSD_db_nofog_visu;

    PSD_db_nofog_nocue_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit_db(1, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_nocue_ap_fit  = vertcat(PSD_db_nofog_nocue_ap_fit{:});
    PSD_db_nofog_audi_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit_db(2, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_audi_ap_fit  = vertcat(PSD_db_nofog_audi_ap_fit{:});
    PSD_db_nofog_visu_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit_db(3, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_visu_ap_fit  = vertcat(PSD_db_nofog_visu_ap_fit{:});

    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit(:,:,1) = PSD_db_nofog_nocue_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit(:,:,2) = PSD_db_nofog_audi_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit(:,:,3) = PSD_db_nofog_visu_ap_fit;

    %%% fog %%%
    PSD_db_fog = PSDtmp(fog_idx);
    PSD_db_fog_nocue_walking     = vertcat(arrayfun(@(s) s.spectra_fooof_db(1, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_nocue_walking     = vertcat(PSD_db_fog_nocue_walking{:});
    PSD_db_fog_nocue_standing    = vertcat(arrayfun(@(s) s.spectra_standing_fooof_db(1, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_nocue_standing    = vertcat(PSD_db_fog_nocue_standing{:});
    PSD_db_fog_nocue             = PSD_db_fog_nocue_walking - PSD_db_fog_nocue_standing;

    PSD_db_fog_audi_walking  = vertcat(arrayfun(@(s) s.spectra_fooof_db(2, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_audi_walking  = vertcat(PSD_db_fog_audi_walking{:});
    PSD_db_fog_audi_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof_db(2, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_audi_standing = vertcat(PSD_db_fog_audi_standing{:});
    PSD_db_fog_audi          = PSD_db_fog_audi_walking - PSD_db_fog_audi_standing;

    PSD_db_fog_visu_walking  = vertcat(arrayfun(@(s) s.spectra_fooof_db(3, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_visu_walking  = vertcat(PSD_db_fog_visu_walking{:});
    PSD_db_fog_visu_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof_db(3, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_visu_standing = vertcat(PSD_db_fog_visu_standing{:});
    PSD_db_fog_visu          = PSD_db_fog_visu_walking - PSD_db_fog_visu_standing;

    PSD_organized.(clusterLabel{clu}).subjects_fog = {PSDtmp(fog_idx).subject}';
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,1) = PSD_db_fog_nocue;
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,2) = PSD_db_fog_audi;
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,3) = PSD_db_fog_visu;

    PSD_db_fog_nocue_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit_db(1, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_nocue_ap_fit  = vertcat(PSD_db_fog_nocue_ap_fit{:});
    PSD_db_fog_audi_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit_db(2, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_audi_ap_fit  = vertcat(PSD_db_fog_audi_ap_fit{:});
    PSD_db_fog_visu_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit_db(3, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_visu_ap_fit  = vertcat(PSD_db_fog_visu_ap_fit{:});

    PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit(:,:,1) = PSD_db_fog_nocue_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit(:,:,2) = PSD_db_fog_audi_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit(:,:,3) = PSD_db_fog_visu_ap_fit;

    % periodic power in uV^2/Hz, walk minus stand
    freqs_uv2 = PSDtmp(1).freqs_fooof;
    groupStructs = {PSD_db_hc, PSD_db_nofog, PSD_db_fog};
    groupFields  = {'PSD_uv2_hc','PSD_uv2_nofog','PSD_uv2_fog'};
    for gi = 1:3
        for cue = 1:3
            tmp = arrayfun(@(s) s.spectra_fooof_power(cue,:) - s.spectra_standing_fooof_power(cue,:), ...
                           groupStructs{gi}, 'UniformOutput', false);
            PSD_organized.(clusterLabel{clu}).(groupFields{gi})(:,:,cue) = vertcat(tmp{:});
        end
    end

    % remove subjects missing ANY cue condition (a missing standing run makes the
    % walk-minus-stand contrast NaN for that cue, not just for visual)
    % compute the removal masks BEFORE any deletion, so the same mask is
    % applied to the data, ap_fit, and subject arrays (they stay aligned)
    rm_hc    = any(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,1,:)),3);
    rm_nofog = any(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,1,:)),3);
    rm_fog   = any(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,1,:)),3);

    disp(['MISSING Cue Conditions for ' clusterLabel{clu}]);
    PSD_organized.(clusterLabel{clu}).subjects_hc(rm_hc)
    PSD_organized.(clusterLabel{clu}).subjects_nofog(rm_nofog)
    PSD_organized.(clusterLabel{clu}).subjects_fog(rm_fog)

    PSD_organized.(clusterLabel{clu}).PSD_db_hc(rm_hc,:,:) = [];
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(rm_nofog,:,:) = [];
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(rm_fog,:,:) = [];

    PSD_organized.(clusterLabel{clu}).PSD_uv2_hc(rm_hc,:,:) = [];
    PSD_organized.(clusterLabel{clu}).PSD_uv2_nofog(rm_nofog,:,:) = [];
    PSD_organized.(clusterLabel{clu}).PSD_uv2_fog(rm_fog,:,:) = [];

    PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit(rm_hc,:,:) = [];
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit(rm_nofog,:,:) = [];
    PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit(rm_fog,:,:) = [];

    PSD_organized.(clusterLabel{clu}).subjects_hc(rm_hc) = [];
    PSD_organized.(clusterLabel{clu}).subjects_nofog(rm_nofog) = [];
    PSD_organized.(clusterLabel{clu}).subjects_fog(rm_fog) = [];

    % merge conditions
    PSD_organized.(clusterLabel{clu}).PSD_db_hc_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_hc,3,'omitmissing');
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_nofog,3,'omitmissing');
    PSD_organized.(clusterLabel{clu}).PSD_db_fog_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_fog,3,'omitmissing');

    PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit,3,'omitmissing');
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit,3,'omitmissing');
    PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit,3,'omitmissing');

end

%% stats - cluster-based permutation

% cue-contrast definitions for the interaction difference scores
% each row: {label, [levelA levelB]}  ->  diff = level A - level B
% dim3 order is 1=nocue, 2=auditory, 3=visual
cueContrasts = { 'audi_nocue', [2 1]; ...
                 'visu_nocue', [3 1]; ...
                 'visu_audi' , [3 2] };

for clu = 1:length(clusterLabel)

    hc    = PSD_organized.(clusterLabel{clu}).PSD_db_hc;      % [nsubj x nfreq x 3cue]
    nofog = PSD_organized.(clusterLabel{clu}).PSD_db_nofog;
    fog   = PSD_organized.(clusterLabel{clu}).PSD_db_fog;
    nfreq = size(hc,2);

    %% ---- Group MAIN effect (omnibus, cue-averaged) : indepsamplesF ----
    hc_coll    = squeeze(mean(hc, 3));       % [nsubj x nfreq]
    nofog_coll = squeeze(mean(nofog, 3));
    fog_coll   = squeeze(mean(fog, 3));

    [gF, ~, gP] = statcondfieldtrip({hc_coll'; nofog_coll'; fog_coll'}, ...
        'paired','off', 'method','permutation', ...
        'naccu',naccu, 'alpha',alpha, 'mcorrect','cluster');
    gF = gF(:); gP = gP(:);
    PSD_organized.(clusterLabel{clu}).groupEffect_pval        = gP;
    [PSD_organized.(clusterLabel{clu}).groupEffect_fval_masked, ~] = apply_significance_mask(gP, gF);

    %% ---- Group PAIRWISE (cue-averaged) : indepsamplesT ----
    % 1=hc vs nofog, 2=hc vs fog, 3=nofog vs fog
    group_pval        = nan(nfreq,3);
    group_stat        = nan(nfreq,3);
    group_tval_masked = nan(nfreq,3);
    groupPairs = {hc_coll nofog_coll; hc_coll fog_coll; nofog_coll fog_coll};
    for p = 1:3
        [group_stat(:,p), ~, group_pval(:,p)] = statcondfieldtrip( ...
            {groupPairs{p,1}'; groupPairs{p,2}'}, ...
            'paired','off', 'method','permutation', ...
            'naccu',naccu, 'alpha',alpha, 'mcorrect','cluster');
        group_pval(:,p) = min(2*group_pval(:,p), 1);
        [group_tval_masked(:,p), ~] = apply_significance_mask(group_pval(:,p), group_stat(:,p), alpha);
    end
    PSD_organized.(clusterLabel{clu}).GroupComp_pval = group_pval;
    PSD_organized.(clusterLabel{clu}).GroupComp_tval = group_tval_masked;

    %% ---- Cue MAIN effect (omnibus, group-pooled) : depsamplesFunivariate ----
    nocue_coll = [hc(:,:,1); nofog(:,:,1); fog(:,:,1)];   % [nsubj_all x nfreq]
    audi_coll  = [hc(:,:,2); nofog(:,:,2); fog(:,:,2)];
    visu_coll  = [hc(:,:,3); nofog(:,:,3); fog(:,:,3)];

    [cF, ~, cP] = statcondfieldtrip({nocue_coll'; audi_coll'; visu_coll'}, ...
        'paired','on', 'method','permutation', ...
        'naccu',naccu, 'alpha',alpha, 'mcorrect','cluster');
    cF = cF(:); cP = cP(:);
    PSD_organized.(clusterLabel{clu}).cueEffect_pval        = cP;
    [PSD_organized.(clusterLabel{clu}).cueEffect_fval_masked, ~] = apply_significance_mask(cP, cF);

    %% ---- Cue PAIRWISE (group-pooled) : depsamplesT ----
    % 1=nocue vs audi, 2=nocue vs visu, 3=audi vs visu
    cue_pval        = nan(nfreq,3);
    cue_stat        = nan(nfreq,3);
    cue_tval_masked = nan(nfreq,3);
    cuePairs = {nocue_coll audi_coll; nocue_coll visu_coll; audi_coll visu_coll};
    for p = 1:3
        [cue_stat(:,p), ~, cue_pval(:,p)] = statcondfieldtrip( ...
            {cuePairs{p,1}'; cuePairs{p,2}'}, ...
            'paired','on', 'method','permutation', ...
            'naccu',naccu, 'alpha',alpha, 'mcorrect','cluster');
        cue_pval(:,p) = min(2*cue_pval(:,p), 1);
        [cue_tval_masked(:,p), ~] = apply_significance_mask(cue_pval(:,p), cue_stat(:,p), alpha);
    end
    PSD_organized.(clusterLabel{clu}).CueComp_pval = cue_pval;
    PSD_organized.(clusterLabel{clu}).CueComp_tval = cue_tval_masked;

    %% ---- Group x cue INTERACTION via difference scores ----
    % For each cue contrast, build per-subject difference spectra by group,
    % then test whether that (within-subject) cue effect differs BETWEEN
    % groups. Omnibus = one-way indepsamplesF across 3 groups; follow-up =
    % pairwise indepsamplesT. All cluster-corrected over frequency.
    nContrast = size(cueContrasts,1);
    inter_omni_pval        = nan(nfreq, nContrast);
    inter_omni_fval_masked = nan(nfreq, nContrast);
    inter_pair_pval        = nan(nfreq, 3, nContrast);   % dim2: group pair
    inter_pair_tval_masked = nan(nfreq, 3, nContrast);
    % store the group-mean difference spectra for plotting
    inter_diff = struct();

    for c = 1:nContrast
        lv = cueContrasts{c,2};                 % [A B] -> A - B
        d_hc    = hc(:,:,lv(1))    - hc(:,:,lv(2));      % [nsubj x nfreq]
        d_nofog = nofog(:,:,lv(1)) - nofog(:,:,lv(2));
        d_fog   = fog(:,:,lv(1))   - fog(:,:,lv(2));

        inter_diff.(cueContrasts{c,1}).hc    = d_hc;
        inter_diff.(cueContrasts{c,1}).nofog = d_nofog;
        inter_diff.(cueContrasts{c,1}).fog   = d_fog;

        % omnibus interaction: does the cue-contrast diff differ across groups?
        [iF, ~, iP] = statcondfieldtrip({d_hc'; d_nofog'; d_fog'}, ...
            'paired','off', 'method','permutation', ...
            'naccu',naccu, 'alpha',alpha, 'mcorrect','cluster');
        iF = iF(:); iP = iP(:);
        inter_omni_pval(:,c) = iP;
        [inter_omni_fval_masked(:,c), ~] = apply_significance_mask(iP, iF);

        % pairwise interaction follow-ups (between-groups on the diff score)
        dPairs = {d_hc d_nofog; d_hc d_fog; d_nofog d_fog};
        for p = 1:3
            [pt, ~, pp] = statcondfieldtrip({dPairs{p,1}'; dPairs{p,2}'}, ...
                'paired','off', 'method','permutation', ...
                'naccu',naccu, 'alpha',alpha, 'mcorrect','cluster');
            inter_pair_pval(:,p,c) = min(2*pp(:), 1);
            [inter_pair_tval_masked(:,p,c), ~] = apply_significance_mask(inter_pair_pval(:,p,c), pt(:), alpha);
        end
    end
    PSD_organized.(clusterLabel{clu}).Interaction_contrasts   = cueContrasts(:,1);
    PSD_organized.(clusterLabel{clu}).Interaction_omni_pval   = inter_omni_pval;
    PSD_organized.(clusterLabel{clu}).Interaction_omni_fval   = inter_omni_fval_masked;
    PSD_organized.(clusterLabel{clu}).Interaction_pair_pval   = inter_pair_pval;
    PSD_organized.(clusterLabel{clu}).Interaction_pair_tval   = inter_pair_tval_masked;
    PSD_organized.(clusterLabel{clu}).Interaction_diff        = inter_diff;

end

%% Group plot  (identical layout to A12; main-effect marker now cluster-based)
%
% Significance markers, drawn as NaN-masked marker rows above the spectra:
%   *  at y = 1.0  omnibus effect (group main effect here, group x cue in figure 2)
%   d  at y = 1.4  post-hoc HC vs PD-NF
%   s  at y = 1.8  post-hoc HC vs PD-F
%   ^  at y = 2.2  post-hoc PD-NF vs PD-F
% Post-hoc rows are gated on the omnibus: if the omnibus has no significant cluster
% anywhere in a panel, that panel's post-hoc rows are blanked, so a pairwise effect
% is never shown without the omnibus that licenses reading it.

hc_color    = [0.15, 0.95, 0.55]; % green
nofog_color = [0.45, 0.00, 0.70]; % purple
fog_color   = [0.00, 0.45, 0.70]; % blue

postHocMarkers = {'d','s','^'};
postHocY       = [1.4 1.8 2.2];
postHocLabels  = {'HC vs PD-NF','HC vs PD-F','PD-NF vs PD-F'};

freqs_plot = PSDtmp(1).freqs_fooof;

f1 = figure('units','normalized','outerposition',[0 0 .99 .33]);
set(gcf, 'Color', 'w','defaultLegendAutoUpdate','off');
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
clusterLabel = {'LSM' 'RSM'};

ax1 = gobjects(2,1);
for clu = 1:length(clusterLabel)
    S = PSD_organized.(clusterLabel{clu});

    ax1(clu) = nexttile(); hold on
    xlim([4 50])
    yline(0,'--')
    xlabel('Frequency (Hz)')
    if clu == 1, ylabel({'Spectral Power', 'Walk - Stand (dB)'}); end

    y = S.PSD_db_hc_combined;
    p1 = plot(freqs_plot,mean(y),'Color',hc_color,'LineWidth',2.6,'LineStyle',':');
    plot_ts_mean_ci(y,freqs_plot,hc_color,[]);

    y = S.PSD_db_nofog_combined;
    p2 = plot(freqs_plot,mean(y),'Color',nofog_color,'LineWidth',2.6,'LineStyle','--');
    plot_ts_mean_ci(y,freqs_plot,nofog_color,[]);

    y = S.PSD_db_fog_combined;
    p3 = plot(freqs_plot,mean(y),'Color',fog_color,'LineWidth',2.7,'LineStyle','-');
    plot_ts_mean_ci(y,freqs_plot,fog_color,[]);

    [pMain, mainDrawn] = plot_sig_marker(freqs_plot, S.groupEffect_pval, alpha, 1, '*', true);
    showPostHoc = ~gatePostHocOnOmnibus || any(S.groupEffect_pval < alpha);
    pPost = gobjects(1,3);
    postDrawn = false(1,3);
    for p = 1:3
        [pPost(p), postDrawn(p)] = plot_sig_marker(freqs_plot, S.GroupComp_pval(:,p), alpha, ...
            postHocY(p), postHocMarkers{p}, showPostHoc);
    end

    hc_n    = ['HC (n=' num2str(size(S.PSD_db_hc_combined,1)) ')'];
    nofog_n = ['PD-NF (n=' num2str(size(S.PSD_db_nofog_combined,1)) ')'];
    fog_n   = ['PD-F (n=' num2str(size(S.PSD_db_fog_combined,1)) ')'];

    % key only the symbols this panel actually drew
    mainLabel  = {'group main effect'};
    sigHandles = [pMain(mainDrawn) pPost(postDrawn)];
    sigLabels  = [mainLabel(mainDrawn) postHocLabels(postDrawn)];
    if isempty(sigHandles)
        lgdFont = 19; lgdCols = 1;
    else
        lgdFont = 14; lgdCols = 2;
    end

    legend(ax1(clu), [p1 p2 p3 sigHandles], ...
        [{hc_n nofog_n fog_n} sigLabels], ...
        'FontSize',lgdFont, 'NumColumns',lgdCols, 'Location','southeast')
end

linkaxes(ax1)
set(ax1, 'FontSize', 16)
%sgtitle('Group Differences (Cue Condition Averaged) by Cluster','fontsize',22)

% save
filename1 = 'Figure-PSD-cluster.png';
saveas(f1,fullfile(figFolder,filename1));

%% Interaction plot  (difference-score spectra by group + omnibus interaction mask)
% rows = cue contrasts (audi-nocue, visu-nocue, visu-audi), cols = LSM/RSM
contrastTitles = {'Auditory - No cue', 'Visual - No cue', 'Visual - Auditory'};
nContrast = size(cueContrasts,1);

f2 = figure('units','normalized','outerposition',[0 0 .7 .99]);
set(gcf, 'Color', 'w','defaultLegendAutoUpdate','off');
tiledlayout(nContrast,2, 'TileSpacing','compact','Padding','compact');

ax2 = gobjects(nContrast*2,1);
% one legend serves all panels, so track which symbols are drawn anywhere in the figure
sigDrawn2   = false(nContrast*2,4);
sigHandles2 = gobjects(nContrast*2,4);
groupHandles2 = gobjects(1,3);
k = 0;
for c = 1:nContrast
    for clu = 1:length(clusterLabel)
        k = k+1;
        ax2(k) = nexttile(); hold on
        xlim([4 50]); yline(0,'--')
        dd = PSD_organized.(clusterLabel{clu}).Interaction_diff.(cueContrasts{c,1});

        y = dd.hc;
        q1 = plot(freqs_plot,mean(y),'Color',hc_color,'LineWidth',2.4,'LineStyle',':');
        plot_ts_mean_ci(y,freqs_plot,hc_color,[]);
        y = dd.nofog;
        q2 = plot(freqs_plot,mean(y),'Color',nofog_color,'LineWidth',2.4,'LineStyle','--');
        plot_ts_mean_ci(y,freqs_plot,nofog_color,[]);
        y = dd.fog;
        q3 = plot(freqs_plot,mean(y),'Color',fog_color,'LineWidth',2.5,'LineStyle','-');
        plot_ts_mean_ci(y,freqs_plot,fog_color,[]);

        % omnibus interaction (group x this cue contrast), then post-hocs gated on it
        omniP = PSD_organized.(clusterLabel{clu}).Interaction_omni_pval(:,c);
        [qMain, mainDrawn] = plot_sig_marker(freqs_plot, omniP, alpha, 1, '*', true);

        pP = PSD_organized.(clusterLabel{clu}).Interaction_pair_pval(:,:,c);
        showPostHoc = ~gatePostHocOnOmnibus || any(omniP < alpha);
        qPost = gobjects(1,3);
        postDrawn = false(1,3);
        for p = 1:3
            [qPost(p), postDrawn(p)] = plot_sig_marker(freqs_plot, pP(:,p), alpha, ...
                postHocY(p), postHocMarkers{p}, showPostHoc);
        end
        sigHandles2(k,:) = [qMain qPost];
        sigDrawn2(k,:)   = [mainDrawn postDrawn];

        if c == 1, title(clusterLabel{clu}); end
        if clu == 1, ylabel({contrastTitles{c}, 'Cue diff (dB)'}); end
        if c == nContrast, xlabel('Frequency (Hz)'); end
        if k == 2, groupHandles2 = [q1 q2 q3]; end
    end
end

% key only the symbols drawn somewhere in the figure, using tile 2's handles
anyDrawn2 = any(sigDrawn2,1);
sigLabels2 = [{'interaction omnibus'} postHocLabels];
legend(ax2(2), [groupHandles2 sigHandles2(2,anyDrawn2)], ...
    [{'HC','PD-NF','PD-F'} sigLabels2(anyDrawn2)], ...
    'FontSize',11, 'NumColumns',1+any(anyDrawn2), 'Location','southeast')

linkaxes(ax2)
set(ax2,'FontSize',12)

filename2 = 'Figure-PSD-interaction.png';
saveas(f2,fullfile(figFolder,filename2));

%% reporting
disp('######')
disp('REPORT')
disp('######')
disp(' ')

if gatePostHocOnOmnibus
    notPlotted = ' [omnibus n.s. - not plotted]';
else
    notPlotted = ' [omnibus n.s.]';
end

groupPairLabels = {'HC vs PD-NF', 'HC vs PD-F', 'PD-NF vs PD-F'};
cuePairLabels   = {'no-cue vs auditory', 'no-cue vs visual', 'auditory vs visual'};
freqs = PSDtmp(1).freqs_fooof;

for clu = 1:length(clusterLabel)
    L = clusterLabel{clu};
    S = PSD_organized.(L);

    % ---- main group effect (omnibus, cluster-corrected) ----
    report_clusters(sprintf('%s main effect of GROUP (indepsamplesF)', L), ...
        S.groupEffect_fval_masked, S.groupEffect_pval, freqs, 'F');

    % ---- group pairwise post-hoc (cluster-corrected) ----
    groupNote = '';
    if ~any(S.groupEffect_pval < alpha), groupNote = notPlotted; end
    for p = 1:3
        report_clusters(sprintf('  post-hoc GROUP pairwise (%s)', groupPairLabels{p}), ...
            S.GroupComp_tval(:,p), S.GroupComp_pval(:,p), freqs, '|T|', groupNote);
    end
    fprintf('\n');

    % ---- main cue effect (omnibus, cluster-corrected) ----
    report_clusters(sprintf('%s main effect of CUE (depsamplesFunivariate)', L), ...
        S.cueEffect_fval_masked, S.cueEffect_pval, freqs, 'F');

    % ---- cue pairwise post-hoc (cluster-corrected) ----
    % not plotted anywhere, so the tag only flags the null omnibus
    cueNote = '';
    if ~any(S.cueEffect_pval < alpha), cueNote = ' [omnibus n.s.]'; end
    for p = 1:3
        report_clusters(sprintf('  post-hoc CUE pairwise (%s)', cuePairLabels{p}), ...
            S.CueComp_tval(:,p), S.CueComp_pval(:,p), freqs, '|T|', cueNote);
    end
    fprintf('\n');

    % ---- group x cue interaction (per contrast, omnibus + pairwise) ----
    for c = 1:nContrast
        report_clusters(sprintf('%s GROUP x CUE interaction (%s, indepsamplesF on diff scores)', ...
            L, contrastTitles{c}), ...
            S.Interaction_omni_fval(:,c), S.Interaction_omni_pval(:,c), freqs, 'F');
        interNote = '';
        if ~any(S.Interaction_omni_pval(:,c) < alpha), interNote = notPlotted; end
        for p = 1:3
            report_clusters(sprintf('  post-hoc interaction pairwise (%s, %s)', ...
                contrastTitles{c}, groupPairLabels{p}), ...
                S.Interaction_pair_tval(:,p,c), S.Interaction_pair_pval(:,p,c), freqs, '|T|', interNote);
        end
    end
    fprintf('\n');
end

%% local functions

function report_clusters(prefix, maskedStat, pvals, freqs, statLabel, note)
% draws one row of significance markers above the spectra, taking the frequency axis (freqs), a p-value vector (pvals), 
% the significance threshold (alpha), the row's height (ypos),  the marker glyph (marker), and a gate flag (showMarkers, false blanks the row),
% and returning the line handle h so the row can appear in the figure's symbol legend even when nothing is significant.

if nargin < 6, note = ''; end

maskedStat = maskedStat(:);
pvals      = pvals(:);
iv = logical2interval(~isnan(maskedStat)');

if isempty(iv)
    fprintf('%s: no significant cluster.\n', prefix);
    return
end

for k = 1:size(iv,1)
    bins = iv(k,1):iv(k,2);
    fprintf('%s: significant cluster %d/%d, %.1f-%.1f Hz, peak %s = %.2f, p_cluster = %.3f.%s\n', ...
        prefix, k, size(iv,1), freqs(bins(1)), freqs(bins(end)), ...
        statLabel, max(abs(maskedStat(bins))), max(pvals(bins)), note);
end
end

function [h, drawn] = plot_sig_marker(freqs, pvals, alpha, ypos, marker, showMarkers)
% One NaN-masked marker row at height ypos; showMarkers=false blanks the row.
% drawn says whether any marker landed, so the caller can omit it from the legend.

y = nan(numel(pvals),1);
if showMarkers
    y(pvals(:) < alpha) = ypos;
end
drawn = any(~isnan(y));
h = plot(freqs, y, marker, 'MarkerSize', 8, 'Color', [0 0 0], 'LineStyle', 'none');
end
