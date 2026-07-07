%% Targeted Clustering and ERSP Selection
%
% This script performs a targeted clustering approach to identify Independent Components (ICs)
% from EEG dipole fitting results that fall within predefined ROI. Targets: Left and Right Sensorimotor areas
%
%   1. Loads subject ERSP data (already baseline-normalized in A6) from a specified folder (e.g., '06_ERSP_*').
%   2. Includes ICs whose centrotype dipole falls inside the Left/Right sensorimotor ROI ellipsoid (MNI).
%   3. Removes ICs whose scalp topography is a group-level outlier (|r| < 0.3 to a subject-balanced group reference map).
%   4. When more than one IC remains in an ROI for a subject, selects the one with the strongest
%      contra-specific (contralateral minus ipsilateral swing) beta suppression.
%   5. Saves the selected IC's ERSP per subject with its dipole and subject info.
%   6. Outputs a summary table of selected ICs and their distances from ROIs.
%
%
% LM 071426
%%
clear all;
if ispc
    mainDir    = 'C:\clab\DoD-Gait';
elseif ismac
    mainDir    = '/Users/Leo/Git/DoD-Gait';
else
    mainDir    = '/home/leo/Documents/DoD-Gait';
end

% load in libraries
funcPath    = fullfile(mainDir,'code','01_eeg_processing_pipeline','func');
addpath(funcPath)
eeglab; ft_defaults; clc; close all

% set up paths
dataDir     = fullfile(mainDir,'data');
figDir      = fullfile(mainDir,'reports','data_quality');

% prompt user for folder
disp('Load in 06_ERSP folder.')
workingDir = uigetdir(dataDir,'Load in 06_ERSP folder.');

[~, folderName] = fileparts(workingDir);
label       = extractAfter(workingDir,'ERSP_');
saveFolder = fullfile(mainDir, 'data', ['07_selected', repmat(['_' label], 1, ~isempty(label))]);

if ~exist(saveFolder, 'dir')
    mkdir(saveFolder)
else
    prompt = "Delete previous folder? Y/N [N]: ";
    txt = input(prompt,"s");
    if strcmpi(txt,'Y')
        rmdir(saveFolder, 's'); % Delete the folder
        disp('Overwrite previous folder');
        mkdir(saveFolder)
    end
end
subFiles = dir(fullfile(workingDir,'sub*'));

% load demo data
T_demo = readtable(fullfile(dataDir,"demo_data_analysis.xlsx"));

% set up MNI corrdinates for ROI
LSM_roi = [-30 -20 58];
RSM_roi = [30 -20 58];
thres = 30; % radius threshold for cluster

roiInclusion = 'centroid';  % which dipole point must fall inside the ROI ellipsoid:
                            %   'centroid' = centrotype point only
                            %   'median'   = per-run cloud median only (ignores centrotype)
                            %   'either'   = centrotype OR cloud median
roiInclusion = validatestring(roiInclusion, {'centroid','median','either'});
roiInside = @(inCentro,inMedian) (strcmp(roiInclusion,'centroid') && inCentro) || ...
                                 (strcmp(roiInclusion,'median')   && inMedian) || ...
                                 (strcmp(roiInclusion,'either')   && (inCentro || inMedian));

% experimental ellipsoid
radii_thres = [30, 25, 30];  % elongated medial-laterally instead of a sphere
ellipsoidMat = diag(radii_thres .^ 2);  % 3x3 matrix for the ellipsoid

selectedFigDir = fullfile(figDir,'A7_selected_ersp');
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

    % update dipfit coord and rv if there is cluster dipfit data
    if isfield(dipfit,'cltr')
        disp('cluster dipfit detected...using cluster data.')
        for i = 1:length(dipfit.model)
            dipfit.model(i).posxyz = dipfit.cltr(i).posxyz;
            dipfit.model(i).rv = dipfit.cltr(i).rv;
        end
        icawinv = cat(2,dipfit.cltr(:).icawinv);
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

    %%% repersentiveIC-to-ROI vector (kept for distance reporting) %%%
    dist_list_LSM       = xyz - LSM_roi;
    dist_list_RSM       = xyz - RSM_roi;

    %%% point-based ROI overlap: include an IC when the point(s) selected %%%
    % by roiInclusion (centrotype and/or per-run cloud median) fall inside
    % the ROI ellipsoid.
    includeInLSM  = false(length(dipfit.model),1);
    includeInRSM  = false(length(dipfit.model),1);
    median_xyz    = nan(length(dipfit.model),3);        % median of the per-run cloud, per IC
    for c = 1:length(dipfit.model)
        if isfield(dipfit,'cltr') && ~isempty(dipfit.cltr(c).pos_runs)
            all_runs = dipfit.cltr(c).pos_runs;            % [nRun x 3] per-run dipole cloud
        else
            all_runs = xyz(c,:);                           % fallback: centrotype point
        end
        all_runs = all_runs(all(~isnan(all_runs),2),:);          % drop runs without this IC

        repXYZ   = xyz(c,:);                     % centrotype point
        medianXYZ       = median(all_runs, 1);             % median of the per-run cloud
        median_xyz(c,:) = medianXYZ;

        % inside-ellipsoid test for the two representative points; which
        % point(s) count toward inclusion is set by roiInclusion above
        centrotypeOffsetLSM = repXYZ - LSM_roi;
        medianOffsetLSM     = medianXYZ - LSM_roi;
        inLSM_centro = sum((centrotypeOffsetLSM / ellipsoidMat) .* centrotypeOffsetLSM) <= 1;
        inLSM_median = sum((medianOffsetLSM / ellipsoidMat) .* medianOffsetLSM) <= 1;
        includeInLSM(c) = roiInside(inLSM_centro, inLSM_median);

        centrotypeOffsetRSM = repXYZ - RSM_roi;
        medianOffsetRSM     = medianXYZ - RSM_roi;
        inRSM_centro = sum((centrotypeOffsetRSM / ellipsoidMat) .* centrotypeOffsetRSM) <= 1;
        inRSM_median = sum((medianOffsetRSM / ellipsoidMat) .* medianOffsetRSM) <= 1;
        includeInRSM(c) = roiInside(inRSM_centro, inRSM_median);
    end
    min_idx_list_LSM    = find(includeInLSM); % centrotype (and optionally median) inside the ROI
    min_idx_list_RSM    = find(includeInRSM);

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
        subjectRows(end+1, :) = {subject, 'LSM', min_idx_list_LSM(i), norm(dist_list_LSM(min_idx_list_LSM(i),:)), norm(median_xyz(min_idx_list_LSM(i),:) - LSM_roi), xyz(min_idx_list_LSM(i),:), median_xyz(min_idx_list_LSM(i),:)};

        counter = length(LSM) + 1;

        LSM(counter).subject        = subject;
        LSM(counter).index          = min_idx_list_LSM(i);
        LSM(counter).icawinv        = icawinv(:,min_idx_list_LSM(i));
        LSM(counter).chanlocs       = chanlocs;
        LSM(counter).times          = times;
        LSM(counter).freqs          = freqs;
        LSM(counter).warpVals       = warpVals;
        LSM(counter).dist_from_roi  = norm(dist_list_LSM(min_idx_list_LSM(i),:));
        % extract ersp
        tmp = tfdata(min_idx_list_LSM(i)).ersp_selfBase_pChange_mean;
        tmpCat = cat(3, tmp{:});
        LSM(counter).ersp_selfBase_pChange_mean = mean(tmpCat, 3, "omitmissing");

        LSM_dip(counter).posxyz     =  dipfit.model(min_idx_list_LSM(i)).posxyz; % for plotting
        LSM_dip(counter).momxyz     =  dipfit.model(min_idx_list_LSM(i)).momxyz;
        LSM_dip(counter).rv         =  dipfit.model(min_idx_list_LSM(i)).rv;

    end

    for i = 1:length(min_idx_list_RSM)
        subjectRows(end+1, :) = {subject , 'RSM', min_idx_list_RSM(i), norm(dist_list_RSM(min_idx_list_RSM(i),:)), norm(median_xyz(min_idx_list_RSM(i),:) - RSM_roi), xyz(min_idx_list_RSM(i),:), median_xyz(min_idx_list_RSM(i),:)};

        counter = length(RSM) + 1;
        RSM(counter).subject = subject;
        RSM(counter).index = min_idx_list_RSM(i);
        RSM(counter).icawinv = icawinv(:,min_idx_list_RSM(i));
        RSM(counter).chanlocs = chanlocs;
        RSM(counter).dist_from_roi  = norm(dist_list_RSM(min_idx_list_RSM(i),:));
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
finalTable = cell2table(allRows, 'VariableNames', {'Subject', 'Cluster', 'DipoleIndex', 'DistanceFromROI', 'ClusterMedianDistanceFromROI', 'CentrotypeXYZ', 'MedianXYZ'});
LSM_rows = finalTable(finalTable.Cluster == "LSM", :);
RSM_rows = finalTable(finalTable.Cluster == "RSM", :);

%% filter by scalp projection
r_thres = 0.3;
nBoot = [];

[LSM_group_map, LSM_topo_corr, LSM_icawinv_norm, LSM_topo_subject] = getTopoCorr(LSM,nBoot,false);
[RSM_group_map, RSM_topo_corr, RSM_icawinv_norm, RSM_topo_subject] = getTopoCorr(RSM,nBoot,false);

% inspect dist of the r values
f0 = figure('Name','A7 topo r distribution','Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; histogram(LSM_topo_corr,20); hold on; 
xline(-1*r_thres,'r','LineWidth',1.5);
xline(r_thres,'r','LineWidth',1.5);
title(sprintf('LSM (thr=%.3f)',r_thres)); 
xlabel('topo r'); 
ylabel('count');
nexttile; histogram(RSM_topo_corr,20); hold on; 
xline(-1*r_thres,'r','LineWidth',1.5);
xline(r_thres,'r','LineWidth',1.5);
title(sprintf('RSM (thr=%.3f)',r_thres)); 
xlabel('topo r');
exportgraphics(f0, fullfile(figDir,'A7_topo_r_distribution.png'), 'Resolution', 300);

% quality check: for each ROI show the scalp map whose |topo r| lands closest
% to r = 0.1, 0.2, ... 0.9
r_targets = 0.1:0.1:0.9;
f_qc = figure('units','normalized','outerposition',[0 0 1 1]);
tiledlayout(2,numel(r_targets),'TileSpacing','compact','Padding','compact');
for t = 1:numel(r_targets)
    [~, idx] = min(abs(LSM_topo_corr - r_targets(t)));
    nexttile(t);
    topoplot(LSM_icawinv_norm(:,idx), LSM(idx).chanlocs);
    title(sprintf('LSM ~%.1f (r=%.2f)', r_targets(t), LSM_topo_corr(idx)));
end
for t = 1:numel(r_targets)
    [~, idx] = min(abs(RSM_topo_corr - r_targets(t)));
    nexttile(numel(r_targets)+t);
    topoplot(RSM_icawinv_norm(:,idx), RSM(idx).chanlocs);
    title(sprintf('RSM ~%.1f (r=%.2f)', r_targets(t), RSM_topo_corr(idx)));
end
sgtitle('Topo map by correlation level','FontSize',24);
set(gcf,'Color','w');
exportgraphics(f_qc, fullfile(figDir,'A7_topo_r_levels.png'), 'Resolution', 300);
close all

LSM_outliers = abs(LSM_topo_corr) < r_thres;
LSM_outliers_idx = find(LSM_outliers);
disp('LSM topo outliers (subject, r):')
disp([LSM_topo_subject(LSM_outliers_idx)', num2cell(LSM_topo_corr(LSM_outliers_idx)')])

RSM_outliers = abs(RSM_topo_corr) < r_thres;
RSM_outliers_idx = find(RSM_outliers);
disp('RSM topo outliers (subject, r):')
disp([RSM_topo_subject(RSM_outliers_idx)', num2cell(RSM_topo_corr(RSM_outliers_idx)')])

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
centroid1 = median(LSM_icawinv_norm, 2);
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
centroid1 = median(RSM_icawinv_norm, 2);
centroid2 = median(RSM_icawinv_norm(:, RSM_outliers), 2);
centroid3 = median(RSM_icawinv_norm(:, ~RSM_outliers), 2);
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
ylim([4 50]); yticks(f(find(ismember(f, [4 8 13 30 50]))))
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
ylim([4 50]); yticks(f(find(ismember(f, [4 8 13 30 50]))))
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
ylim([4 50]); yticks(f(find(ismember(f, [4 8 13 30 50]))))
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
ylim([4 50]); yticks(f(find(ismember(f, [4 8 13 30 50]))))
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
ylim([4 50]); yticks(f(find(ismember(f, [4 8 13 30 50]))))
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
ylim([4 50]); yticks(f(find(ismember(f, [4 8 13 30 50]))))
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
subFiles = dir(fullfile(workingDir,'sub*'));
delete(gcp('nocreate'))
parpool("Processes",7)
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
        T(sub).LSM_cluster_median_distance_from_roi = clusterMedianDist(ERSP.dipfit,idx_select,LSM_roi);

        if isfield(ERSP,'removedTrialNum')
            T(sub).LSM_cycles_removed           = ERSP.removedTrialNum{idx_select};
            T(sub).LSM_cycles_removed_percentage= ERSP.removedTrialNum{idx_select}/size(ERSP.data,4)*100;
        end

        saveName = fullfile(saveFolder,['LSM_' ERSP.subject '.mat']);
        parsave(saveName,ERSP,idx_select,'LSM')
    end

    %%% get RSM IC %%%
    if ~isempty(min_idx_list_RSM)
        disp(['RSM IC found for subject: ' ERSP.subject])

        idx_select = selectIC(ERSP,min_idx_list_RSM,'RSM',selectedFigDir);

        T(sub).RSM_ic_number                = idx_select;
        T(sub).RSM_distance_from_roi        = pdist2(ERSP.dipfit.model(idx_select).posxyz,RSM_roi);
        T(sub).RSM_cluster_median_distance_from_roi = clusterMedianDist(ERSP.dipfit,idx_select,RSM_roi);

        if isfield(ERSP,'removedTrialNum')
            T(sub).RSM_cycles_removed           = ERSP.removedTrialNum{idx_select};
            T(sub).RSM_cycles_removed_percentage= ERSP.removedTrialNum{idx_select}/size(ERSP.data,4)*100;
        end

        saveName = fullfile(saveFolder,['RSM_' ERSP.subject '.mat']);
        parsave(saveName,ERSP,idx_select,'RSM')
    end
end
delete(gcp('nocreate'))

%% save to file
tableName = fullfile(figDir, 'A7_cluster_data_table.csv');
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
    % get the ic with the strongest contra-specific beta suppression
    % (contralateral minus ipsilateral swing), so a lateralized sensorimotor source
    contra_swing    = squeeze(mean(data(:,f_bound(1):f_bound(2),contra_startIdx:contra_endIdx),2));
    ips_swing       = squeeze(mean(data(:,f_bound(1):f_bound(2),ips_startIdx:ips_endIdx),2));
    contra_beta     = mean(contra_swing,2);
    ips_beta        = mean(ips_swing,2);
    lat_suppression = contra_beta - ips_beta; % more negative = stronger contra suppression
    
    [~, idx] = min(lat_suppression);

    % contra_swing    = squeeze(mean(data(:,f_bound(1):f_bound(2),contra_startIdx:contra_endIdx),2));
    % contra_beta_suppression = sum(contra_swing,2);
    % [~, idx] = min(contra_beta_suppression); % get the greatest beta suppression
    %
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
%%% if the IC is in both clusters, assign it to the hemisphere contralateral to the swing with stronger beta suppression

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

function d = clusterMedianDist(dipfit,ic,roi)
%%% Euclidean distance from the cluster's median per-run dipole position to
% an ROI. Falls back to the centrotype position when no per-run cloud exists.
if isfield(dipfit,'cltr') && ~isempty(dipfit.cltr(ic).pos_runs)
    runs = dipfit.cltr(ic).pos_runs;
    runs = runs(all(~isnan(runs),2),:);   % drop runs without this IC
    medianXYZ = median(runs,1);
else
    medianXYZ = dipfit.model(ic).posxyz;
end
d = pdist2(medianXYZ,roi);
end