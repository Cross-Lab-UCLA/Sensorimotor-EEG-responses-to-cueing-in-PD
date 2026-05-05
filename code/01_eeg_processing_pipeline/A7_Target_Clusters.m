%% Targeted Clustering and ERSP Normalization
%
% This script performs a targeted clustering approach to identify Independent Components (ICs)
% from EEG dipole fitting results that are closest to predefined Regions of
% Interest (ROIs). Targets: Left and Right Sensorimotor areas
%
%   1. Loads subject ERSP data from a specified folder (e.g., '06_ERSP_subjects_*').
%   2. Identifies IC dipoles nearest to ROIs defined by MNI coordinates for Left and Right Sensorimotor areas.
%   3. Resolves overlapping ICs found in both ROIs by keeping the closest one.
%   4. Normalizes ERSP data of selected ICs with respect to baseline and across experimental conditions.
%   5. Saves normalized ERSP data and clustering information for each subject.
%   6. Outputs a summary table of selected ICs and their distances from ROIs.
%
% LM 122325
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
eeglab; ft_defaults; clc; close all

% set up paths
dataDir     = fullfile(mainDir,'data');
figDir      = fullfile(mainDir,'reports','data_quality');

% prompt user for folder
disp('Load in 06_ERSP folder.')
% dataDir2     = 'C:\clab\DoD-Gait\data';
path = uigetdir(dataDir,'Load in 06_ERSP folder.');
label = extractAfter(path,'06_ERSP_');
saveDir     = fullfile(mainDir,'data',['07_selected_' label]);
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
else
    prompt = "Delete previous folder? Y/N [N]: ";
    txt = input(prompt,"s");
    if strcmpi(txt,'Y')
        rmdir(saveDir, 's'); % Delete the folder
        disp('Overwrite previous folder');
        mkdir(saveDir)
    end
end
subFiles = dir(fullfile(path,'sub*'));

% load demo data
T_demo = readtable(fullfile(dataDir,"demo_data_analysis.xlsx"));

% set up MNI corrdinates for ROI
LSM_roi = [-30 -20 58];
RSM_roi = [30 -20 58];
thres = 30; % radius threshold for cluster

% experimental ellipsoid
radii_thres = [30, 25, 30];  % elongated medial-laterally instead of a sphere
ellipsoidMat = diag(radii_thres .^ 2);  % 3x3 matrix for the ellipsoid

selectedFigDir = fullfile(figDir,['selected_ersp_' label]);
if ~exist(selectedFigDir, 'dir')
    mkdir(selectedFigDir)
end

%% Check clusters
allSubjectsSelectedDipole = cell(length(subFiles), 1);
LSM = [];
RSM = [];
LSM_icawinv_subject_mean = nan(64,length(subFiles));
RSM_icawinv_subject_mean = nan(64,length(subFiles));

for sub = 1:length(subFiles)
    disp(['Loading ' subFiles(sub).name])
    load(fullfile(subFiles(sub).folder,subFiles(sub).name));
    disp('Checking dipole... ')

    % update dipfit corrd and rv if there is bootstrapped data
    if isfield(dipfit,'boot')
        disp('bootstrap dipfit detected...using bootrapped data.')
        for i = 1:length(dipfit.model)
            dipfit.model(i).posxyz = dipfit.boot(i).posxyz;
            dipfit.model(i).momxyz = dipfit.boot(i).momxyz;
            dipfit.model(i).rv = dipfit.boot(i).rv;
            icawinv = cat(2,dipfit.boot(:).icawinv);
        end
    end

    % Make sure channel order are consistent
    if sub == 1
        labels1 = {chanlocs.labels};
    else
        labels2 = {chanlocs.labels};
        [~, idx] = ismember(labels1, labels2);
        chanlocs = chanlocs(idx);
        icawinv = icawinv(idx,:);
    end
    
    % get corrdinates of dipoles
    xyz = nan(length(dipfit.model),3);
    for c = 1:size(xyz,1) % transfer posiiton and rv data to local variables
        xyz(c,:) = dipfit.model(c).posxyz;
    end

    %%% find distance to LSM and RSM using sphere %%%
    % dist_list_LSM = pdist2(xyz,LSM_roi);
    % min_idx_list_LSM = find(dist_list_LSM < thres);
    % dist_list_RSM = pdist2(xyz,RSM_roi);
    % min_idx_list_RSM = find(dist_list_RSM < thres);

    %%% find distance to LSM and RSM using ellipsoid %%%
    dist_list_LSM       = xyz - LSM_roi;
    min_idx_list_LSM    = find(sum((dist_list_LSM / ellipsoidMat) .* dist_list_LSM, 2) <= 1); % find the dipoles within ellipsoid
    dist_list_RSM       = xyz - RSM_roi;
    min_idx_list_RSM    = find(sum((dist_list_RSM / ellipsoidMat) .* dist_list_RSM, 2) <= 1); % find the dipoles within ellipsoid

    % find if there are shared IC and remove it from the cluster with the
    % furthest dist
    % shared_ic = intersect(min_idx_list_LSM, min_idx_list_RSM);
    % if ~isempty(shared_ic)
    %     disp('Shared IC found...selecting according to contralateral beta...')
    %     ERSP = load(fullfile(subFiles(sub).folder,subFiles(sub).name));
    %     [min_idx_list_LSM, min_idx_list_RSM] = checkShared(ERSP,shared_ic,min_idx_list_LSM,min_idx_list_RSM);
    %     ERSP = [];
    % end

    subjectRows = {};
    for i = 1:length(min_idx_list_LSM)
        subjectRows(end+1, :) = {subject, 'LSM', min_idx_list_LSM(i), dist_list_LSM(min_idx_list_LSM(i)), xyz(min_idx_list_LSM(i),:)};

        counter = length(LSM) + 1;

        LSM(counter).subject        = subject;
        LSM(counter).index          = min_idx_list_LSM(i);
        LSM(counter).icawinv        = icawinv(:,min_idx_list_LSM(i));
        LSM(counter).chanlocs       = chanlocs;
        LSM(counter).times          = times;
        LSM(counter).freqs          = freqs;
        LSM(counter).warpVals       = warpVals;
        LSM(counter).dist_from_roi  = dist_list_LSM(min_idx_list_LSM(i));
        % extract ersp
        tmp = tfdata(min_idx_list_LSM(i)).ersp_selfBase_pChange_mean;
        tmpCat = cat(3, tmp{:});
        LSM(counter).ersp_selfBase_pChange_mean = mean(tmpCat, 3, "omitmissing");

        LSM_dip(counter).posxyz     =  dipfit.model(min_idx_list_LSM(i)).posxyz;
        LSM_dip(counter).momxyz     =  dipfit.model(min_idx_list_LSM(i)).momxyz;
        LSM_dip(counter).rv         =  dipfit.model(min_idx_list_LSM(i)).rv;

    end

    for i = 1:length(min_idx_list_RSM)
        subjectRows(end+1, :) = {subject , 'RSM', min_idx_list_RSM(i), dist_list_RSM(min_idx_list_RSM(i)), xyz(min_idx_list_RSM(i),:)};

        counter = length(RSM) + 1;
        RSM(counter).subject = subject;
        RSM(counter).index = min_idx_list_RSM(i);
        RSM(counter).icawinv = icawinv(:,min_idx_list_RSM(i));
        RSM(counter).chanlocs = chanlocs;
        RSM(counter).dist_from_roi  = dist_list_RSM(min_idx_list_RSM(i));
        RSM(counter).times          = times;
        RSM(counter).freqs          = freqs;
        RSM(counter).warpVals       = warpVals;
        % extract ersps
        tmp = tfdata(min_idx_list_RSM(i)).ersp_selfBase_pChange_mean;
        tmpCat = cat(3, tmp{:});
        RSM(counter).ersp_selfBase_pChange_mean = mean(tmpCat, 3,"omitmissing");

        RSM_dip(counter).posxyz   =  dipfit.model(min_idx_list_RSM(i)).posxyz;
        RSM_dip(counter).momxyz   =  dipfit.model(min_idx_list_RSM(i)).momxyz;
        RSM_dip(counter).rv   =  dipfit.model(min_idx_list_RSM(i)).rv;
    end

    allSubjectsSelectedDipole{sub} = subjectRows;

end
allRows = vertcat(allSubjectsSelectedDipole{:});
finalTable = cell2table(allRows, 'VariableNames', {'Subject', 'Cluster', 'DipoleIndex', 'DistanceFromROI', 'XYZ'});
LSM_rows = finalTable(finalTable.Cluster == "LSM", :);
RSM_rows = finalTable(finalTable.Cluster == "RSM", :);

%% filter by scalp projection
r_thres = 0.3; % dipoles with r value below this is removed
nBoot = [];

[LSM_group_map LSM_topo_corr,LSM_icawinv_norm] = getTopoCorr(LSM,nBoot);
LSM_outliers = LSM_topo_corr < r_thres;
LSM_outliers_idx = find(LSM_outliers);

[RSM_group_map RSM_topo_corr,RSM_icawinv_norm] = getTopoCorr(RSM,nBoot);
RSM_outliers = RSM_topo_corr < r_thres;
RSM_outliers_idx = find(RSM_outliers);

% LSM plot outliers
f1 = figure('units','normalized','outerposition',[0 0 1 1]);
all_vals = LSM_icawinv_norm(:, LSM_outliers_idx);
clim = [min(all_vals(:)), max(all_vals(:))];
for out = 1:length(LSM_outliers_idx)
    nexttile;
    topoplot(LSM_icawinv_norm(:, LSM_outliers_idx(out)), LSM(LSM_outliers_idx(out)).chanlocs, ...
        'maplimits', clim);
    title([LSM(LSM_outliers_idx(out)).subject])
end
set(gcf, 'Color', 'w');
figName = fullfile(figDir, 'A7_Removed_Dipoles_LSM.png');
exportgraphics(gcf, figName, 'Resolution', 500);
close all

% RSM plot outliers
f1 = figure('units','normalized','outerposition',[0 0 1 1]);
all_vals = RSM_icawinv_norm(:, RSM_outliers_idx);
clim = [min(all_vals(:)), max(all_vals(:))];
for out = 1:length(RSM_outliers_idx)
    nexttile;
    topoplot(RSM_icawinv_norm(:, RSM_outliers_idx(out)), RSM(RSM_outliers_idx(out)).chanlocs, ...
        'maplimits', clim);
    title([RSM(RSM_outliers_idx(out)).subject])
end
figName = fullfile(figDir, 'A7_Removed_Dipoles_RSM.png');
exportgraphics(gcf, figName, 'Resolution', 500);
close all

%% plot before and after topo filtering %%
f1 = figure('units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(2,4,'TileSpacing', 'compact', 'Padding', 'compact');

% the get the channel locs for plotting average (plotting purposes only)
% these subject has a more uniformed chan pos on the topo
LSM_plot_idx = find(ismember({LSM.subject},{'sub-HC15','sub-HC17','sub-PD12','sub-PD23','sub-PD31'}));
RSM_plot_idx = find(ismember({RSM.subject},{'sub-HC15','sub-HC17','sub-PD12','sub-PD23','sub-PD31'}));

% LSM
centroid1 = LSM_group_map;
centroid2 = median(LSM_icawinv_norm(:, LSM_outliers), 2);
centroid3 = median(LSM_icawinv_norm(:, ~LSM_outliers), 2);
all_vals1 = [centroid1(:); centroid2(:); centroid3(:)];
clim1 = [min(all_vals1), max(all_vals1)];

LSM_num = size(LSM_icawinv_norm,2);
LSM_num_removed = length(find(LSM_outliers));

nexttile(tl,1);
topoplot(centroid1,LSM(LSM_plot_idx(1)).chanlocs, 'maplimits', clim1);
%title(['LSM Before n=' num2str(LSM_num)]);

nexttile(tl,2);
try
    topoplot(centroid2,LSM(LSM_plot_idx(1)).chanlocs, 'maplimits', clim1);
    %title(['LSM Outliers n=' num2str(LSM_num_removed)]);
end

nexttile(tl,3);
try
    topoplot(centroid3,LSM(LSM_plot_idx(1)).chanlocs, 'maplimits', clim1);
    %title(['LSM After n=' num2str(LSM_num - LSM_num_removed)]);
end

nexttile(tl,4); axis off;  % blank tile for colorbar
cb1 = colorbar;
cb1.Layout.Tile = 4;
cb1.Label.String = 'LSM Scalp Map';

% RSM
centroid1 = RSM_group_map;
centroid2 = mean(RSM_icawinv_norm(:, RSM_outliers), 2);
centroid3 = mean(RSM_icawinv_norm(:, ~RSM_outliers), 2);
all_vals2 = [centroid1(:); centroid2(:); centroid3(:)];
clim2 = [min(all_vals2), max(all_vals2)];

RSM_num = size(RSM_icawinv_norm,2);
RSM_num_removed = length(find(RSM_outliers));

nexttile(tl,5);
topoplot(centroid1,RSM(RSM_plot_idx(1)).chanlocs, 'maplimits', clim2);
%title(['RSM Before n=' num2str(RSM_num)]);

nexttile(tl,6);
try
    topoplot(centroid2,RSM(RSM_plot_idx(1)).chanlocs, 'maplimits', clim2);
    %title(['RSM Outliers n=' num2str(RSM_num_removed)]);
end

nexttile(tl,7);
try
    topoplot(centroid3,RSM(RSM_plot_idx(1)).chanlocs, 'maplimits', clim2);
    %title(['RSM After n=' num2str(RSM_num - RSM_num_removed)]);
end

nexttile(tl,8); axis off;  % blank tile for colorbar
cb2 = colorbar;
cb2.Layout.Tile = 8;
cb2.Label.String = 'RSM Scalp Map';

sgtitle('Scalp Map Filtering','FontSize',36);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 32);
set(gcf, 'Color', 'w');

figName = fullfile(figDir, 'A7_Scalp_Map_Cluster_Cleaning.png');
exportgraphics(gcf, figName, 'Resolution', 500);
close all

%% plot by ERSP before and after
LSM_ersp_vec = cat(3, LSM(:).ersp_selfBase_pChange_mean);
RSM_ersp_vec = cat(3, RSM(:).ersp_selfBase_pChange_mean);

f1 = figure('units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(2,10,'TileSpacing', 'compact', 'Padding', 'compact');

ersp_og = mean(LSM_ersp_vec, 3, 'omitnan');
in_idx  = setdiff(1:size(LSM_ersp_vec,3), LSM_outliers_idx);
ersp_in = mean(LSM_ersp_vec(:,:, in_idx), 3, 'omitnan');
ersp_out = mean(LSM_ersp_vec(:,:, LSM_outliers_idx), 3, 'omitnan');
t = LSM(1).times./1000;
f = LSM(1).freqs;
warpVals = LSM(1).warpVals./1000;

nexttile(tl,1,[1 3]);
p1 = contourf(t,f,ersp_og,40,'linecolor','none');hold on
set(gca, 'ydir', 'normal'); ylabel('Frequencies (Hz)','FontSize',18);
xlim([0 warpVals(end)]); xticks(warpVals);
xline(warpVals([1,3,5]),'LineWidth',2,'alpha',1);
xline(warpVals([2,4]),'LineWidth',2,'alpha',1);
ylim([4 50]); yticks(freqs(find(ismember(freqs, [4 8 13 30 50]))))
xticklabels({'RHS','LTO','LHS','RTO','RHS'});xtickangle(30);
xlabel('Gait Cycle')
yline(13,'--k',LineWidth=2)
yline(30,'--k',LineWidth=2)
title(['LSM mean ERSP before outlier removal n=' num2str(size(LSM_ersp_vec,3))]);

nexttile(tl,4,[1 3]);
p2 = contourf(t,f,ersp_out,40,'linecolor','none');hold on
set(gca, 'ydir', 'normal'); ylabel('Frequencies (Hz)','FontSize',18);
xlim([0 warpVals(end)]); xticks(warpVals);
xline(warpVals([1,3,5]),'LineWidth',2,'alpha',1);
xline(warpVals([2,4]),'LineWidth',2,'alpha',1);
ylim([4 50]); yticks(freqs(find(ismember(freqs, [4 8 13 30 55]))))
xticklabels({'RHS','LTO','LHS','RTO','RHS'});xtickangle(30);
xlabel('Gait Cycle')
yline(13,'--k',LineWidth=2)
yline(30,'--k',LineWidth=2)
title(['LSM mean outlier ERSP removed n=' num2str(length(LSM_outliers_idx))]);

nexttile(tl,7,[1 3]);
p3 = contourf(t,f,ersp_in,40,'linecolor','none');hold on
set(gca, 'ydir', 'normal'); ylabel('Frequencies (Hz)','FontSize',18);
xlim([0 warpVals(end)]); xticks(warpVals);
xline(warpVals([1,3,5]),'LineWidth',2,'alpha',1);
xline(warpVals([2,4]),'LineWidth',2,'alpha',1);
ylim([4 50]); yticks(freqs(find(ismember(freqs, [4 8 13 30 55]))))
xticklabels({'RHS','LTO','LHS','RTO','RHS'});xtickangle(30);
xlabel('Gait Cycle')
yline(13,'--k',LineWidth=2)
yline(30,'--k',LineWidth=2)
title(['LSM mean ERSP kept n=' num2str(length(in_idx))]);

nexttile(tl,10); % right side of grid
axis off
cb1 = colorbar;
cb1.Label.String = 'LSM ERSP Color Map';
allVals = [ersp_og(:); ersp_in(:); ersp_out(:)];
climVals = [min(allVals)+1 max(allVals)-1];
ax = findall(f1,'type','axes');
set(ax,'CLim',climVals);

% RSM
ersp_og = mean(RSM_ersp_vec, 3, 'omitnan');
in_idx  = setdiff(1:size(RSM_ersp_vec,3), RSM_outliers_idx);
ersp_in = mean(RSM_ersp_vec(:,:, in_idx), 3, 'omitnan');
ersp_out = mean(RSM_ersp_vec(:,:, RSM_outliers), 3, 'omitnan');

nexttile(tl,11,[1 3]);
p1 = contourf(t,f,ersp_og,40,'linecolor','none');hold on
set(gca, 'ydir', 'normal'); ylabel('Frequencies (Hz)','FontSize',18);
xlim([0 warpVals(end)]); xticks(warpVals);
xline(warpVals([1,3,5]),'LineWidth',2,'alpha',1);
xline(warpVals([2,4]),'LineWidth',2,'alpha',1);
ylim([4 50]); yticks(freqs(find(ismember(freqs, [4 8 13 30 55]))))
xticklabels({'RHS','LTO','LHS','RTO','RHS'});xtickangle(30);
xlabel('Gait Cycle')
yline(13,'--k',LineWidth=2)
yline(30,'--k',LineWidth=2)
title(['RSM mean ERSP before outlier removal n=' num2str(size(RSM_ersp_vec,3))]);

nexttile(tl,14,[1 3]);
p2 = contourf(t,f,ersp_out,40,'linecolor','none');hold on
set(gca, 'ydir', 'normal'); ylabel('Frequencies (Hz)','FontSize',18);
xlim([0 warpVals(end)]); xticks(warpVals);
xline(warpVals([1,3,5]),'LineWidth',2,'alpha',1);
xline(warpVals([2,4]),'LineWidth',2,'alpha',1);
ylim([4 50]); yticks(freqs(find(ismember(freqs, [4 8 13 30 55]))))
xticklabels({'RHS','LTO','LHS','RTO','RHS'});xtickangle(30);
xlabel('Gait Cycle')
yline(13,'--k',LineWidth=2)
yline(30,'--k',LineWidth=2)
title(['RSM mean outlier ERSP removed n=' num2str(length(RSM_outliers_idx))]);

nexttile(tl,17,[1 3]);
p3 = contourf(t,f,ersp_in,40,'linecolor','none');hold on
set(gca, 'ydir', 'normal'); ylabel('Frequencies (Hz)','FontSize',18);
xlim([0 warpVals(end)]); xticks(warpVals);
xline(warpVals([1,3,5]),'LineWidth',2,'alpha',1);
xline(warpVals([2,4]),'LineWidth',2,'alpha',1);
ylim([4 50]); yticks(freqs(find(ismember(freqs, [4 8 13 30 55]))))
xticklabels({'RHS','LTO','LHS','RTO','RHS'});xtickangle(30);
xlabel('Gait Cycle')
yline(13,'--k',LineWidth=2)
yline(30,'--k',LineWidth=2)
title(['RSM mean ERSP kept n=' num2str(length(in_idx))]);

nexttile(tl,20);
axis off
cb1 = colorbar;
cb1.Label.String = 'RSM ERSP Color Map';
allVals = [ersp_og(:); ersp_in(:); ersp_out(:)];
climVals = [min(allVals)+1 max(allVals)-1];
ax = findall(f1,'type','axes');
set(ax,'CLim',climVals);
sgtitle('ERSP Filtering','FontSize',30);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 16);
set(gcf, 'Color', 'w');
figName = fullfile(figDir, 'A7_Scalp_Map_Cluster_Cleaning_ERSP.png');
exportgraphics(gcf, figName, 'Resolution', 500);
%close all

%% remove outliers from topo matching %%
LSM(LSM_outliers) = [];
RSM(RSM_outliers) = [];

%% plot dipole
colorRSM1 = {[.9 .43 .3]};
colorLSM1 = {[.3 .82 .9]};
colorOut = {[.9 .9 .9]};

f1 = figure('units','normalized','outerposition',[0 0 1 1]);
tl = tiledlayout(2,4,'TileSpacing', 'compact', 'Padding', 'compact');
views = [45 45; 0 90; 90 0; 0 0];

LSM_dip_kept = LSM_dip;
LSM_dip_kept(LSM_outliers) = [];
LSM_dip_all = [LSM_dip LSM_dip_kept];
for x = 1:4
    nexttile
    colorList = [repmat(colorLSM1, length(LSM_dip_all), 1)];
    colorList(1:length(LSM_dip)) = colorOut;
    dipplot(LSM_dip_all, ...
        'gui','off',...
        'verbose','off', ...
        'spheres','on', ...
        'dipolelength',0,...
        'view', views(x,:), ...
        'color',colorList,...
        'dipolesize', 25*ones(length(LSM_dip_all),1));

    axis vis3d equal
    camzoom(0.85)         % Keep zoom consistent
    set(gca, 'XLim', [-100 100], 'YLim', [-100 100], 'ZLim', [-100 100]);  % consistent space
end

RSM_dip_kept = RSM_dip;
RSM_dip_kept(RSM_outliers) = [];
RSM_dip_all = [RSM_dip RSM_dip_kept];

for x = 1:4
    nexttile
    colorList = [repmat(colorRSM1, length(RSM_dip_all), 1)];
    colorList(1:length(RSM_dip)) = colorOut;
    dipplot(RSM_dip_all, ...
        'gui','off',...
        'verbose','off', ...
        'spheres','on', ...
        'dipolelength',0,...
        'view', views(x,:), ...
        'color',colorList,...
        'dipolesize', 25*ones(length(RSM_dip_all),1));

    axis vis3d equal
    camzoom(0.85)         % Keep zoom consistent
    set(gca, 'XLim', [-100 100], 'YLim', [-100 100], 'ZLim', [-100 100]);  % consistent space
end

figName = fullfile(figDir, 'A7_Dipole_Cluster_AfterCleaning.png');
exportgraphics(gcf, figName, 'Resolution', 500);
close all

%% run sorting loop
subFiles = dir(fullfile(path,'sub*'));
delete(gcp('nocreate'))
parpool("Processes",8)
parfor sub = 1:length(subFiles)

    % load EEG
    disp(['Selecting ERSP: ' subFiles(sub).name]);
    ERSP = load(fullfile(subFiles(sub).folder,subFiles(sub).name));

    % add subgroup
    ERSP.subgroup   = T_demo.subgroup(strcmp(T_demo.subject,ERSP.subject));
    ERSP.height     = T_demo.Height(strcmp(T_demo.subject,ERSP.subject));
    ERSP.weight     = T_demo.Weight(strcmp(T_demo.subject,ERSP.subject));
    ERSP.age        = T_demo.Age(strcmp(T_demo.subject,ERSP.subject));
    ERSP.MOCA       = T_demo.MOCA(strcmp(T_demo.subject,ERSP.subject));
    ERSP.UPDRS3     = T_demo.UPDRS3(strcmp(T_demo.subject,ERSP.subject));
    ERSP.PAS        = T_demo.PAS(strcmp(T_demo.subject,ERSP.subject));

    % add info to table
    T(sub).subject  = ERSP.subject;
    T(sub).group    = ERSP.group;
    T(sub).gaitCycleUsed = length(ERSP.HS_idx.kept4TF) / ...
        (length(ERSP.HS_idx.kept4TF) + length(ERSP.HS_idx.removed4TF));

    LSM_idx = find(strcmp({LSM.subject},ERSP.subject));
    RSM_idx = find(strcmp({RSM.subject},ERSP.subject));

    %%% get LSM and RSM list %%%
    min_idx_list_LSM  = [LSM(LSM_idx).index];
    min_idx_list_RSM  = [RSM(RSM_idx).index];

    shared_ic = intersect(min_idx_list_LSM, min_idx_list_RSM);
    if ~isempty(shared_ic)
        disp('Shared IC found...selecting according to contralateral beta...')
        [min_idx_list_LSM, min_idx_list_RSM] = checkShared(ERSP,shared_ic,min_idx_list_LSM,min_idx_list_RSM);
    end

    %%% example plot
    % figure; hold on;
    % [u, v] = meshgrid(linspace(0, 2*pi, 50), linspace(0, pi, 50));
    % scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 20, 'k');  % All dipoles
    % scatter3(xyz(min_idx_list_RSM,1), xyz(min_idx_list_RSM,2), xyz(min_idx_list_RSM,3), 40, 'r', 'filled');  % LSM dipoles
    % scatter3(xyz(min_idx_list_LSM,1), xyz(min_idx_list_LSM,2), xyz(min_idx_list_LSM,3), 40, 'b', 'filled');  % RSM dipoles
    % R_x = radii_thres(1) * cos(u) .* sin(v) + RSM_roi(1);
    % R_y = radii_thres(2) * sin(u) .* sin(v) + RSM_roi(2);
    % R_z = radii_thres(3) * cos(v) + RSM_roi(3);
    % L_x = radii_thres(1) * cos(u) .* sin(v) + LSM_roi(1);
    % L_y = radii_thres(2) * sin(u) .* sin(v) + LSM_roi(2);
    % L_z = radii_thres(3) * cos(v) + LSM_roi(3);
    % surf(R_x, R_y, R_z, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'r');
    % surf(L_x, L_y, L_z, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'b');
    % xlabel('X'); ylabel('Y'); zlabel('Z');
    % title('Dipoles within Left/Right Sensorimotor Ellipsoids');
    % axis equal;
    % grid on;

    %%% get LSM IC %%%
    if ~isempty(min_idx_list_LSM)
        disp(['LSM IC found for subject: ' ERSP.subject])

        idx_select = selectIC(ERSP,min_idx_list_LSM,'LSM',selectedFigDir);

        T(sub).LSM_ic_number            = idx_select;
        T(sub).LSM_distance_from_roi    = pdist2(ERSP.dipfit.model(idx_select).posxyz,LSM_roi);

        if isfield(ERSP,'removedTrialNum')
            T(sub).LSM_cycles_removed           = ERSP.removedTrialNum{idx_select};
            T(sub).LSM_cycles_removed_percentage= ERSP.removedTrialNum{idx_select}/size(ERSP.data,4)*100;
        end

        saveName = fullfile(saveDir,['LSM_' ERSP.subject '.mat']);
        parsave(saveName,ERSP,idx_select,'LSM')
    end

    %%% get RSM IC %%%
    if ~isempty(min_idx_list_RSM)
        disp(['RSM IC found for subject: ' ERSP.subject])

        idx_select = selectIC(ERSP,min_idx_list_RSM,'RSM',selectedFigDir);

        T(sub).RSM_ic_number                = idx_select;
        T(sub).RSM_distance_from_roi        = pdist2(ERSP.dipfit.model(idx_select).posxyz,RSM_roi);

        if isfield(ERSP,'removedTrialNum')
            T(sub).RSM_cycles_removed           = ERSP.removedTrialNum{idx_select};
            T(sub).RSM_cycles_removed_percentage= ERSP.removedTrialNum{idx_select}/size(ERSP.data,4)*100;
        end

        saveName = fullfile(saveDir,['RSM_' ERSP.subject '.mat']);
        parsave(saveName,ERSP,idx_select,'RSM')
    end
end
delete(gcp('nocreate'))

%% save to file
tableName = fullfile(figDir, ['A7_cluster_data_table_' label '.csv']);
if exist(tableName, 'file') == 2
    delete(tableName); % Delete the file
    disp('Overwrite previous csv');
end
writetable(struct2table(T),tableName);
disp('done.')

%% functions %%
%% Function to get the IC with the greatest contralateral beta suppression
% should there be multiple ICs picked up in the roi cluster
function idx_select = selectIC(ERSP,min_idx_list,roi,figDir)

data = [];
for i = 1:length(ERSP.tfdata)
    ersp = ERSP.tfdata(i).ersp_selfBase_pChange_mean;
    ersp = cat(3, ersp{:});
    erspAvg = mean(ersp, 3, 'omitmissing');
    data(:, :, i) = erspAvg;
end
%data = cat(3, ERSP.tfdata(:).ersp_p_allCond_mean);
data = permute(data, [3 1 2]);
data = squeeze(data(min_idx_list,:,:));

switch roi
    case 'LSM'
        contra_startIdx = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(4),'nearest');
        contra_endIdx   = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(5),'nearest');
        ips_startIdx    = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(2),'nearest');
        ips_endIdx      = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(3),'nearest');
        f_bound         = [find(ERSP.freqs == 13) find(ERSP.freqs == 30)];

    case 'RSM'
        contra_startIdx = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(2),'nearest');
        contra_endIdx   = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(3),'nearest');
        ips_startIdx    = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(4),'nearest');
        ips_endIdx      = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(5),'nearest');
        f_bound         = [find(ERSP.freqs == 13) find(ERSP.freqs == 30)];

end

if length(min_idx_list) > 1
    % get the ic with the greatest beta suppression
    contra_swing    = squeeze(mean(data(:,f_bound(1):f_bound(2),contra_startIdx:contra_endIdx),2));
    contra_beta_suppression = sum(contra_swing,2);
    [~, idx] = min(contra_beta_suppression); % get the greatest beta suppression

    % beta_power = squeeze(mean(data(:,f_bound(1):f_bound(2),:),2));
    % total_abs_beta_power = sum(abs(beta_power),2);
    % [~, idx] = min(contra_beta_suppression./total_abs_beta_power); % get the most relative beta suppression

    idx_select = min_idx_list(idx);
    plotERSP(ERSP,data,f_bound,min_idx_list,idx,figDir,roi)

else
    idx_select = min_idx_list;
    plotERSP(ERSP,data,f_bound,min_idx_list,1,figDir,roi)
end
end


function parsave(filename,ERSP,idx_select,roi)
%save function
disp(['Saving ' filename])

% save struct
tfdata = ERSP.tfdata(idx_select);

% add additional info
tfdata.subgroup   = ERSP.subgroup;
tfdata.height     = ERSP.height;
tfdata.weight     = ERSP.weight;
tfdata.age        = ERSP.age;
tfdata.MOCA       = ERSP.MOCA;
tfdata.UPDRS3     = ERSP.UPDRS3;
tfdata.PAS        = ERSP.PAS;
tfdata.roi        = roi;
tfdata.HS_idx     = ERSP.HS_idx;

save(filename,'tfdata','-v7.3','-nocompression')
end

function [min_idx_list_LSM, min_idx_list_RSM] = checkShared(ERSP,shared_ic,min_idx_list_LSM,min_idx_list_RSM)
%%% if the IC is both clusters, keep the IC that is closest to the ROI

for i = 1:length(ERSP.tfdata)
    ersp = ERSP.tfdata(i).ersp_selfBase_pChange_mean;
    ersp = cat(3, ersp{:});
    erspAvg = mean(ersp, 3, 'omitmissing');
    data(:, :, i) = erspAvg;
end
data = permute(data, [3 1 2]);

for sh = 1:length(shared_ic)

    dataS = squeeze(data(shared_ic(sh),:,:));
    L_startIdx = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(2),'nearest');
    L_endIdx   = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(3),'nearest');
    R_startIdx    = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(4),'nearest');
    R_endIdx      = interp1(ERSP.times,1:length(ERSP.times),ERSP.warpVals(5),'nearest');
    f_bound         = [find(ERSP.freqs == 13) find(ERSP.freqs == 30)];

    L_swing    = squeeze(mean(dataS(f_bound(1):f_bound(2),L_startIdx:L_endIdx),1));
    L_swing_beta_suppression = mean(L_swing);

    R_swing    = squeeze(mean(dataS(f_bound(1):f_bound(2),R_startIdx:R_endIdx),1));
    R_swing_beta_suppression = mean(R_swing);

    if R_swing_beta_suppression <= L_swing_beta_suppression
        % if more beta is suppressed during the right swing, then keep
        % the dipole in the LSM, and remove the it from RSM
        min_idx_list_RSM(min_idx_list_RSM ==shared_ic(sh)) = [];
    else
        min_idx_list_LSM(min_idx_list_LSM ==shared_ic(sh)) = [];
    end
end

end