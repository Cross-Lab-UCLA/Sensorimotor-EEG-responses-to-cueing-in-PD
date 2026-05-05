%% Make cued averaged group PSD 
% This script performs PSD analysis and plot
%
% LM 121225
%%
clear all; clc; eeglab; close all
rng(134);
if ispc
    mainDir     = 'E:\clab\DoD-Gait';
elseif ismac
    mainDir = '/Users/Leo/Git/DoD-Gait';
else
    mainDir    = '/home/leo/Documents/DoD-Gait';
end
funcpath    = fullfile(mainDir,'code','func');
addpath(funcpath)

% directories
dataDir = fullfile(mainDir, 'data');
figFolder  = fullfile(mainDir, 'reports', 'results');

% load demo data
T_demo = readtable(fullfile(dataDir,"demo_data_analysis.xlsx"));

% load data
psdFile  = dir(fullfile(dataDir, '08_PSD_withStanding_run1only_moreAffected','*.mat'));
load(fullfile(psdFile.folder,psdFile.name));

%% organize
clusterLabel = {'LSM' 'RSM'};
PSD_organized.LSM = [];
PSD_organized.RSM = [];

for clu = 1:length(clusterLabel)
    idx = [];
    idx = find(strcmp({PSD.roi},clusterLabel{clu}));
    PSDtmp = PSD(idx);
    hc_idx = find(strcmp({PSDtmp.subgroup},'hc'));
    fog_idx = find(strcmp({PSDtmp.subgroup},'fog'));
    nofog_idx = find(strcmp({PSDtmp.subgroup},'nofog'));

    %%% HC %%%
    PSD_db_hc = PSDtmp(hc_idx);
    PSD_db_hc_nocue_walking     = vertcat(arrayfun(@(s) s.spectra_fooof(1, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_nocue_walking     = vertcat(PSD_db_hc_nocue_walking{:});
    PSD_db_hc_nocue_standing    = vertcat(arrayfun(@(s) s.spectra_standing_fooof(1, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_nocue_standing    = vertcat(PSD_db_hc_nocue_standing{:});
    PSD_db_hc_nocue             = PSD_db_hc_nocue_walking - PSD_db_hc_nocue_standing;

    PSD_db_hc_audi_walking  = vertcat(arrayfun(@(s) s.spectra_fooof(2, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_audi_walking  = vertcat(PSD_db_hc_audi_walking{:});
    PSD_db_hc_audi_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof(2, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_audi_standing = vertcat(PSD_db_hc_audi_standing{:});
    PSD_db_hc_audi          = PSD_db_hc_audi_walking - PSD_db_hc_audi_standing;

    PSD_db_hc_visu_walking  = vertcat(arrayfun(@(s) s.spectra_fooof(3, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_visu_walking  = vertcat(PSD_db_hc_visu_walking{:});
    PSD_db_hc_visu_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof(3, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_visu_standing = vertcat(PSD_db_hc_visu_standing{:});
    PSD_db_hc_visu          = PSD_db_hc_visu_walking-PSD_db_hc_visu_standing;

    PSD_organized.(clusterLabel{clu}).subjects_hc = {PSDtmp(hc_idx).subject}';
    PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,1) = PSD_db_hc_nocue;
    PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,2) = PSD_db_hc_audi;
    PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,3) = PSD_db_hc_visu;
    
    PSD_db_hc_nocue_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit(1, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_nocue_ap_fit  = vertcat(PSD_db_hc_nocue_ap_fit{:});
    PSD_db_hc_audi_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit(2, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_audi_ap_fit  = vertcat(PSD_db_hc_audi_ap_fit{:});
    PSD_db_hc_visu_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit(3, :), PSD_db_hc, 'UniformOutput', false));
    PSD_db_hc_visu_ap_fit  = vertcat(PSD_db_hc_visu_ap_fit{:});

    PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit(:,:,1) = PSD_db_hc_nocue_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit(:,:,2) = PSD_db_hc_audi_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit(:,:,3) = PSD_db_hc_visu_ap_fit;

    %%% nofog %%%
    PSD_db_nofog = PSDtmp(nofog_idx);
    PSD_db_nofog_nocue_walking     = vertcat(arrayfun(@(s) s.spectra_fooof(1, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_nocue_walking     = vertcat(PSD_db_nofog_nocue_walking{:});
    PSD_db_nofog_nocue_standing    = vertcat(arrayfun(@(s) s.spectra_standing_fooof(1, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_nocue_standing    = vertcat(PSD_db_nofog_nocue_standing{:});
    PSD_db_nofog_nocue             = PSD_db_nofog_nocue_walking - PSD_db_nofog_nocue_standing;

    PSD_db_nofog_audi_walking  = vertcat(arrayfun(@(s) s.spectra_fooof(2, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_audi_walking  = vertcat(PSD_db_nofog_audi_walking{:});
    PSD_db_nofog_audi_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof(2, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_audi_standing = vertcat(PSD_db_nofog_audi_standing{:});
    PSD_db_nofog_audi          = PSD_db_nofog_audi_walking - PSD_db_nofog_audi_standing;

    PSD_db_nofog_visu_walking  = vertcat(arrayfun(@(s) s.spectra_fooof(3, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_visu_walking  = vertcat(PSD_db_nofog_visu_walking{:});
    PSD_db_nofog_visu_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof(3, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_visu_standing = vertcat(PSD_db_nofog_visu_standing{:});
    PSD_db_nofog_visu          = PSD_db_nofog_visu_walking - PSD_db_nofog_visu_standing;

    PSD_organized.(clusterLabel{clu}).subjects_nofog = {PSDtmp(nofog_idx).subject}';
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,1) = PSD_db_nofog_nocue;
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,2) = PSD_db_nofog_audi;
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,3) = PSD_db_nofog_visu;

    PSD_db_nofog_nocue_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit(1, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_nocue_ap_fit  = vertcat(PSD_db_nofog_nocue_ap_fit{:});
    PSD_db_nofog_audi_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit(2, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_audi_ap_fit  = vertcat(PSD_db_nofog_audi_ap_fit{:});
    PSD_db_nofog_visu_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit(3, :), PSD_db_nofog, 'UniformOutput', false));
    PSD_db_nofog_visu_ap_fit  = vertcat(PSD_db_nofog_visu_ap_fit{:});

    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit(:,:,1) = PSD_db_nofog_nocue_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit(:,:,2) = PSD_db_nofog_audi_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit(:,:,3) = PSD_db_nofog_visu_ap_fit;

    %%% nofog %%%
    PSD_db_fog = PSDtmp(fog_idx);
    PSD_db_fog_nocue_walking     = vertcat(arrayfun(@(s) s.spectra_fooof(1, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_nocue_walking     = vertcat(PSD_db_fog_nocue_walking{:});
    PSD_db_fog_nocue_standing    = vertcat(arrayfun(@(s) s.spectra_standing_fooof(1, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_nocue_standing    = vertcat(PSD_db_fog_nocue_standing{:});
    PSD_db_fog_nocue             = PSD_db_fog_nocue_walking - PSD_db_fog_nocue_standing;

    PSD_db_fog_audi_walking  = vertcat(arrayfun(@(s) s.spectra_fooof(2, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_audi_walking  = vertcat(PSD_db_fog_audi_walking{:});
    PSD_db_fog_audi_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof(2, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_audi_standing = vertcat(PSD_db_fog_audi_standing{:});
    PSD_db_fog_audi          = PSD_db_fog_audi_walking - PSD_db_fog_audi_standing;

    PSD_db_fog_visu_walking  = vertcat(arrayfun(@(s) s.spectra_fooof(3, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_visu_walking  = vertcat(PSD_db_fog_visu_walking{:});
    PSD_db_fog_visu_standing = vertcat(arrayfun(@(s) s.spectra_standing_fooof(3, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_visu_standing = vertcat(PSD_db_fog_visu_standing{:});
    PSD_db_fog_visu          = PSD_db_fog_visu_walking - PSD_db_fog_visu_standing;

    PSD_organized.(clusterLabel{clu}).subjects_fog = {PSDtmp(fog_idx).subject}';
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,1) = PSD_db_fog_nocue;
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,2) = PSD_db_fog_audi;
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,3) = PSD_db_fog_visu;

    PSD_db_fog_nocue_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit(1, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_nocue_ap_fit  = vertcat(PSD_db_fog_nocue_ap_fit{:});
    PSD_db_fog_audi_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit(2, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_audi_ap_fit  = vertcat(PSD_db_fog_audi_ap_fit{:});
    PSD_db_fog_visu_ap_fit  = vertcat(arrayfun(@(s) s.ap_fit(3, :), PSD_db_fog, 'UniformOutput', false));
    PSD_db_fog_visu_ap_fit  = vertcat(PSD_db_fog_visu_ap_fit{:});

    PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit(:,:,1) = PSD_db_fog_nocue_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit(:,:,2) = PSD_db_fog_audi_ap_fit;
    PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit(:,:,3) = PSD_db_fog_visu_ap_fit;

    % remove subjects with no visu

    disp(['MISSING Visual Cue Conditions for ' clusterLabel{clu}]);
    PSD_organized.(clusterLabel{clu}).subjects_hc(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,1,3)))
    PSD_organized.(clusterLabel{clu}).subjects_nofog(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,1,3)))
    PSD_organized.(clusterLabel{clu}).subjects_fog(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,1,3)))

    PSD_organized.(clusterLabel{clu}).PSD_db_hc(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,1,3)),:,:) = [];
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,1,3)),:,:) = [];
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,1,3)),:,:) = [];

    PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,1,3)),:,:) = [];
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,1,3)),:,:) = [];
    PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit(isnan(PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,1,3)),:,:) = [];

    % merge conditions
    PSD_organized.(clusterLabel{clu}).PSD_db_hc_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_hc,3,'omitmissing');
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_nofog,3,'omitmissing');
    PSD_organized.(clusterLabel{clu}).PSD_db_fog_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_fog,3,'omitmissing');

    PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_hc_ap_fit,3,'omitmissing');
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_nofog_ap_fit,3,'omitmissing');
    PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit_combined = mean(PSD_organized.(clusterLabel{clu}).PSD_db_fog_ap_fit,3,'omitmissing');

end

%% stats
for clu = 1:length(clusterLabel)

% [pcond, pgroup, pinter, statcond, statgroup, statinter] = std_stat(...
%     { ...
%         PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,1)'...
%         PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,1)'...
%         PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,1)'; ...
%         ...
%         PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,2)'...
%         PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,2)'...
%         PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,2)'; ...
%         ...
%         PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,3)'...
%         PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,3)'...
%         PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,3)'
%         },...
%     'groupstats','on',...
%     'condstats','on',...
%     'paired',{'on','off'}, ...
%     'fieldtripmethod','montecarlo', ...
%     'fieldtripmcorrect','cluster',...
%     'mode','fieldtrip');

[pcond, pgroup, pinter, stat_cond, stat_group, stat_inter] = std_stat(...
    { ...
        PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,1)'...
        PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,1)'...
        PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,1)'; ...
        ...
        PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,2)'...
        PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,2)'...
        PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,2)'; ...
        ...
        PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,3)'...
        PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,3)'...
        PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,3)'
        },...
    'groupstats','on',...
    'condstats','on',...
    'paired',{'on','off'}, ...
    'method','permutation', ...
    'naccu',2000, ...
    'mcorrect','fdr');

PSD_organized.(clusterLabel{clu}).interaction_pval = pinter{3};
PSD_organized.(clusterLabel{clu}).groupEffect_pval = pinter{1};
PSD_organized.(clusterLabel{clu}).cueEffect_pval = pinter{2};

[inter_masked_fvals, ~] = apply_significance_mask(pinter{3}, stat_inter{3});
[group_masked_fvals, ~] = apply_significance_mask(pinter{1}, stat_inter{1});
[cue_masked_fvals, ~]   = apply_significance_mask(pinter{2}, stat_inter{2});

PSD_organized.(clusterLabel{clu}).interaction_fval_masked  = inter_masked_fvals;
PSD_organized.(clusterLabel{clu}).groupEffect_fval_masked  = group_masked_fvals;
PSD_organized.(clusterLabel{clu}).cueEffect_fval_masked    = cue_masked_fvals;

%%% main group effects with cluster correction
hc_coll   = squeeze(mean(PSD_organized.(clusterLabel{clu}).PSD_db_hc, 3));
nofog_coll = squeeze(mean(PSD_organized.(clusterLabel{clu}).PSD_db_nofog, 3));
fog_coll   = squeeze(mean(PSD_organized.(clusterLabel{clu}).PSD_db_fog, 3));

group_pval = nan(size(hc_coll,2),3);
group_stat = nan(size(hc_coll,2),3);
group_tval_masked = nan(size(hc_coll,2),3);

% 1st col = hc vs nofog
% 2nd col = hc vs fog
% 3rd col = nofog vs fog

[group_stat(:,1),df,group_pval(:,1)]= statcondfieldtrip({hc_coll'; nofog_coll'},...
    'paired','off',...
    'method','permutation',...
    'naccu', 2000,...
    'alpha', 0.05,...
    'mcorrect','cluster');

[group_tval_masked(:,1), ~] = apply_significance_mask(group_pval(:,1), group_stat(:,1))

[group_stat(:,2),df,group_pval(:,2)]= statcondfieldtrip({hc_coll'; fog_coll'},...
    'paired','off',...
    'method','permutation',...
    'naccu', 2000,...
    'alpha', 0.05,...
    'mcorrect','cluster');

[group_tval_masked(:,2), ~] = apply_significance_mask(group_pval(:,2), group_stat(:,2))

[group_stat(:,3),df,group_pval(:,3)]= statcondfieldtrip({nofog_coll'; fog_coll'},...
    'paired','off',...
    'method','permutation',...
    'naccu', 2000,...
    'alpha', 0.05,...
    'mcorrect','cluster');

[group_tval_masked(:,3), ~] = apply_significance_mask(group_pval(:,3), group_stat(:,3))

PSD_organized.(clusterLabel{clu}).GroupComp_pval = group_pval;
PSD_organized.(clusterLabel{clu}).GroupComp_tval = group_tval_masked;

%%% main CUE effects  mcorrect
nocue_coll = [PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,1);...
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,1);...
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,1)];

audi_coll = [PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,2);...
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,2);...
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,2)];

visu_coll = [PSD_organized.(clusterLabel{clu}).PSD_db_hc(:,:,3);...
    PSD_organized.(clusterLabel{clu}).PSD_db_nofog(:,:,3);...
    PSD_organized.(clusterLabel{clu}).PSD_db_fog(:,:,3)];

cue_pval = nan(size(nocue_coll,2),3);
cue_stat = nan(size(hc_coll,2),3);
cue_tval_masked = nan(size(hc_coll,2),3);

% 1st col = nocue vs auditory
% 2nd col = nocue vs visual
% 3rd col = auditory vs visual

[cue_stat(:,1),df,cue_pval(:,1)]= statcondfieldtrip({nocue_coll'; audi_coll'},...
    'paired','on',...
    'method','permutation',...
    'naccu', 2000,...
    'alpha', 0.05,...
    'mcorrect','cluster');
[cue_tval_masked(:,1), ~] = apply_significance_mask(cue_pval(:,1), cue_stat(:,1))

[cue_stat(:,2),df,cue_pval(:,2)]= statcondfieldtrip({nocue_coll'; visu_coll'},...
    'paired','on',...
    'method','permutation',...
    'naccu', 2000,...
    'alpha', 0.05,...
    'mcorrect','cluster');
[cue_tval_masked(:,2), ~] = apply_significance_mask(cue_pval(:,2), cue_stat(:,2))

[cue_stat(:,3),df,cue_pval(:,3)]= statcondfieldtrip({audi_coll'; visu_coll'},...
    'paired','on',...
    'method','permutation',...
    'naccu', 2000,...
    'alpha', 0.05,...
    'mcorrect','cluster');
[cue_tval_masked(:,3), ~] = apply_significance_mask(cue_pval(:,3), cue_stat(:,3))


PSD_organized.(clusterLabel{clu}).CueComp_pval = cue_pval;
PSD_organized.(clusterLabel{clu}).CueComp_tval = cue_tval_masked;

end

%% Group plot

hc_color    = [0.15, 0.95, 0.55]; % green
nofog_color = [0.45, 0.00, 0.70]; % purple
fog_color   = [0.00, 0.45, 0.70]; % blue

f1 = figure('units','normalized','outerposition',[0 0 .99 .33]);
set(gcf, 'Color', 'w','defaultLegendAutoUpdate','off');
main_tl = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
clusterLabel = {'LSM' 'RSM'};

%LSM
ax1(1) = nexttile(); hold on
%title('Left Sensorimotor');
xlim([4 50])
yline(0,'--')
xlabel('Frequency (Hz)')
ylabel({'Spectral Power', 'Walk - Stand (dB)'});

y = PSD_organized.LSM.PSD_db_hc_combined;
p1 = plot(PSDtmp(1).freqs_fooof,mean(y),'Color',hc_color,'LineWidth',2.6,'LineStyle',':'); hold on
plot_ts_mean_ci(y,PSDtmp(1).freqs_fooof,hc_color,[]); hold on

y = PSD_organized.LSM.PSD_db_nofog_combined;
p2 = plot(PSDtmp(1).freqs_fooof,mean(y),'Color',nofog_color,'LineWidth',2.6,'LineStyle','--');
plot_ts_mean_ci(y,PSDtmp(1).freqs_fooof,nofog_color,[]);

y = PSD_organized.LSM.PSD_db_fog_combined;
p3 = plot(PSDtmp(1).freqs_fooof,mean(y),'Color',fog_color,'LineWidth',2.7,'LineStyle','-');
plot_ts_mean_ci(y,PSDtmp(1).freqs_fooof,fog_color,[]);

signif_bar_main = nan(length(PSD_organized.LSM.groupEffect_pval),1);
signif_bar_main(PSD_organized.LSM.groupEffect_pval<0.05) = 1;
plot(PSDtmp(1).freqs_fooof,signif_bar_main*.1,'*','MarkerSize',7,'Color',[0 0 0])

signif_bar_1 = nan(length(PSD_organized.LSM.GroupComp_pval(:,3)),1); % HC vs NOFOG
signif_bar_1(PSD_organized.LSM.GroupComp_pval(:,1)<0.05) = 1;
signif_bar_2 = nan(length(PSD_organized.LSM.GroupComp_pval(:,3)),1); % HC vs FOG
signif_bar_2(PSD_organized.LSM.GroupComp_pval(:,2)<0.05) = 1;
signif_bar_3 = nan(length(PSD_organized.LSM.GroupComp_pval(:,3)),1); % NOFOG vs FOG
signif_bar_3(PSD_organized.LSM.GroupComp_pval(:,3)<0.05) = 1;

sigidx = find(signif_bar_1 == 1);
for i = 1:length(sigidx)
    text(PSDtmp(1).freqs_fooof(sigidx(i)), ... 
         signif_bar_1(sigidx(i)) * 0.14, ...
         '♦', 'HorizontalAlignment','center','FontSize',10);
end

sigidx = find(signif_bar_2 == 1);
for i = 1:length(sigidx)
    text(PSDtmp(1).freqs_fooof(sigidx(i)), ...
         signif_bar_2(sigidx(i)) * 0.18, ...
         'b', 'HorizontalAlignment','center','FontSize',14);
end

sigidx = find(signif_bar_3 == 1);
for i = 1:length(sigidx)
    text(PSDtmp(1).freqs_fooof(sigidx(i)), ...
         signif_bar_3(sigidx(i)) * 0.22, ...
         'c', 'HorizontalAlignment','center','FontSize',14);
end

hc_n = ['HC (n=' num2str(size(PSD_organized.LSM.PSD_db_hc_combined,1)) ')'];
nofog_n = ['PD-NF (n=' num2str(size(PSD_organized.LSM.PSD_db_nofog_combined,1)) ')'];
fog_n = ['PD-F (n=' num2str(size(PSD_organized.LSM.PSD_db_fog_combined,1)) ')'];

legend(ax1(1),[p1 p2 p3], ...
    {hc_n nofog_n fog_n},...
    'FontSize',19,....
    'Location','southeast')

% RSM
ax1(2) = nexttile(); hold on
%title('Right Sensorimotor');
xlim([4 50])
yline(0,'--')
xlabel('Frequency (Hz)')

y = PSD_organized.RSM.PSD_db_hc_combined;
p4 = plot(PSDtmp(1).freqs_fooof,mean(y),'Color',hc_color,'LineWidth',2.6,'LineStyle',':'); hold on
plot_ts_mean_ci(y,PSDtmp(1).freqs_fooof,hc_color,[]); hold on

y = PSD_organized.RSM.PSD_db_nofog_combined;
p5 = plot(PSDtmp(1).freqs_fooof,mean(y),'Color',nofog_color,'LineWidth',2.6,'LineStyle','--');
plot_ts_mean_ci(y,PSDtmp(1).freqs_fooof,nofog_color,[]);

y = PSD_organized.RSM.PSD_db_fog_combined;
p6 = plot(PSDtmp(1).freqs_fooof,mean(y),'Color',fog_color,'LineWidth',2.7,'LineStyle','-');
plot_ts_mean_ci(y,PSDtmp(1).freqs_fooof,fog_color,[]);

signif_bar_main = nan(length(PSD_organized.RSM.groupEffect_pval),1);
signif_bar_main(PSD_organized.RSM.groupEffect_pval<0.05) = 1;
plot(PSDtmp(1).freqs_fooof,signif_bar_main*.1,'*','MarkerSize',7,'Color',[0 0 0])

signif_bar_1 = nan(length(PSD_organized.RSM.GroupComp_pval(:,3)),1); % HC vs NOFOG
signif_bar_1(PSD_organized.RSM.GroupComp_pval(:,1)<0.05) = 1;
signif_bar_2 = nan(length(PSD_organized.RSM.GroupComp_pval(:,3)),1);
signif_bar_2(PSD_organized.RSM.GroupComp_pval(:,2)<0.05) = 1;
signif_bar_3 = nan(length(PSD_organized.RSM.GroupComp_pval(:,3)),1);
signif_bar_3(PSD_organized.RSM.GroupComp_pval(:,3)<0.05) = 1;

sigidx = find(signif_bar_1 == 1);
for i = 1:length(sigidx)

    text(PSDtmp(1).freqs_fooof(sigidx(i)), ... 
         signif_bar_1(sigidx(i)) * 0.14, ...
         '♦', 'HorizontalAlignment','center','FontSize',10);
end

sigidx = find(signif_bar_2 == 1);
for i = 1:length(sigidx)
    text(PSDtmp(1).freqs_fooof(sigidx(i)), ...
         signif_bar_2(sigidx(i)) * 0.18, ...
         'b', 'HorizontalAlignment','center','FontSize',14);
end

sigidx = find(signif_bar_3 == 1);
for i = 1:length(sigidx)
    text(PSDtmp(1).freqs_fooof(sigidx(i)), ...
         signif_bar_3(sigidx(i)) * 0.22, ...
         'c', 'HorizontalAlignment','center','FontSize',14);
end

hc_n = ['HC (n=' num2str(size(PSD_organized.RSM.PSD_db_hc_combined,1)) ')'];
nofog_n = ['PD-NF (n=' num2str(size(PSD_organized.RSM.PSD_db_nofog_combined,1)) ')'];
fog_n = ['PD-F (n=' num2str(size(PSD_organized.RSM.PSD_db_fog_combined,1)) ')'];

legend(ax1(2),[p4 p5 p6], ...
    {hc_n nofog_n fog_n},...
    'FontSize',19,....
    'Location','southeast')


%ylabel('db')
linkaxes(ax1)
set(ax1, 'FontSize', 16)
%sgtitle('Group Differences (Cue Condition Averaged) by Cluster','fontsize',22)

% save
filename1 = 'Figure-PSD.png';
saveas(f1,fullfile(figFolder,filename1));

%% reporting

% sigifcance at main group level for LSM
LSM_F_sum = sum(PSD_organized.LSM.groupEffect_fval_masked,'omitmissing');
LSM_min_F = min(PSD_organized.LSM.groupEffect_fval_masked, [], 'all', 'omitnan');
LSM_max_F = max(PSD_organized.LSM.groupEffect_fval_masked, [], 'all', 'omitnan');
lsm_freqs = PSDtmp(1).freqs_fooof(~isnan(PSD_organized.LSM.groupEffect_fval_masked));
LSM_pvals= PSD_organized.LSM.groupEffect_pval(~isnan(PSD_organized.LSM.groupEffect_fval_masked));

fprintf('For the LSM group, FDR-corrected permutation testing revealed a significant main effect\n');
fprintf('of group in the %.1f–%.1f Hz frequency range.\n', min(lsm_freqs), max(lsm_freqs));
fprintf('Statistics: F-range = [%.2f, %.2f], peak F = %.2f, all p_fdr < %.3f.\n\n', ...
    LSM_min_F, LSM_max_F, LSM_max_F, max(LSM_pvals));

% sigifcance at main group level for RSM
RSM_F_sum = sum(PSD_organized.RSM.groupEffect_fval_masked, 'omitmissing');
RSM_min_F = min(PSD_organized.RSM.groupEffect_fval_masked, [], 'all', 'omitnan');
RSM_max_F = max(PSD_organized.RSM.groupEffect_fval_masked, [], 'all', 'omitnan');
rsm_freqs = PSDtmp(1).freqs_fooof(~isnan(PSD_organized.RSM.groupEffect_fval_masked))
RSM_pvals= PSD_organized.RSM.groupEffect_pval(~isnan(PSD_organized.RSM.groupEffect_fval_masked));

fprintf('For the RSM group, FDR-corrected permutation testing revealed a significant main effect\n');
fprintf('of group in the %.1f–%.1f Hz frequency range.\n', min(rsm_freqs), max(rsm_freqs));
fprintf('Statistics: F-range = [%.2f, %.2f], peak F = %.2f, all p_fdr < %.3f.\n\n', ...
    RSM_min_F, RSM_max_F, RSM_max_F, max(RSM_pvals));


% significance in RSM between hc and nofog
T_sum = sum(PSD_organized.RSM.GroupComp_tval(:,1),'omitmissing')
comp_freqs = PSDtmp(1).freqs_fooof(~isnan(PSD_organized.RSM.GroupComp_tval(:,1)));

fprintf('Post-hoc pairwise comparisons for the RSM group (HC vs. No-FoG) demonstrated\n');
fprintf('a significant cluster in the %.1f–%.1f Hz range (Sum T = %.2f).\n', ...
    min(comp_freqs), max(comp_freqs), T_sum);