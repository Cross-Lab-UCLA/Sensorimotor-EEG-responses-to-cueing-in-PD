%% plot figure of allersp mean-ed
% This script perform group level statistical analysis of the ersp roi 
% and make the ERSP plot for manuscript
%
% LM 122325
%%
clear all; clc; close all
rng(134);
if ispc
    mainDir     = 'E:\clab\DoD-Gait';
else
    mainDir = '/Users/Leo/Git/DOD-Gait';
end
addpath(fullfile(mainDir,'code'));
addpath(fullfile(mainDir,'code','func'));
eeglab;
ft_defaults; close all;

dataFolder = fullfile(mainDir,'data','07_selected_withStanding_run1only_moreAffected');
figFolder = fullfile(mainDir,'reports','results');

%%
rsm_ersp_selfBase_pChange_mean = [];
lsm_ersp_selfBase_pChange_mean = [];
rsm_ersp_allBase_pChange_mean = [];
lsm_ersp_allBase_pChange_mean = [];

clusterLabel = {'RSM' 'LSM'};

rsmFiles = dir(fullfile(dataFolder,'RSM*'));
for s = 1:length(rsmFiles)

    disp(['Loading ' rsmFiles(s).name])
    load(fullfile(rsmFiles(s).folder,rsmFiles(s).name));
    disp('Loading complete.');

    rsm_subject{s}  = tfdata.subject;
    rsm_subgroup(s) = tfdata.subgroup;

    rsm_ersp_mean_db(:,:,s) = tfdata.ersp_db_allCond_mean;
    
    rsm_ersp_selfBase_pChange               = cat(3,tfdata.ersp_selfBase_pChange{:});
    rsm_ersp_selfBase_pChange_mean(:,:,s)   = mean(rsm_ersp_selfBase_pChange,3,"omitnan");

    rsm_ersp_allBase_pChange               = cat(3,tfdata.ersp_allBase_pChange{:});
    rsm_ersp_allBase_pChange_mean(:,:,s)   = mean(rsm_ersp_allBase_pChange,3,"omitnan");

    rsm_dip(s).posxyz   =  tfdata.dipfit.posxyz;
    rsm_dip(s).momxyz   =  tfdata.dipfit.momxyz;
    rsm_dip(s).rv       =  tfdata.dipfit.rv;
end

lsmFiles = dir(fullfile(dataFolder,'LSM*'));
for s = 1:length(lsmFiles)

    disp(['Loading ' lsmFiles(s).name])
    load(fullfile(lsmFiles(s).folder,lsmFiles(s).name));
    disp('Loading complete.');

    lsm_subject{s}  = tfdata.subject;
    lsm_subgroup(s) = tfdata.subgroup;

    lsm_ersp_mean_db(:,:,s) = tfdata.ersp_db_allCond_mean;

    lsm_ersp_selfBase_pChange               = cat(3,tfdata.ersp_selfBase_pChange{:});
    lsm_ersp_selfBase_pChange_mean(:,:,s)   = mean(lsm_ersp_selfBase_pChange,3,"omitnan");

    lsm_ersp_allBase_pChange               = cat(3,tfdata.ersp_allBase_pChange{:});
    lsm_ersp_allBase_pChange_mean(:,:,s)   = mean(lsm_ersp_allBase_pChange,3,"omitnan");

    lsm_dip(s).posxyz   =  tfdata.dipfit.posxyz;
    lsm_dip(s).momxyz   =  tfdata.dipfit.momxyz;
    lsm_dip(s).rv       =  tfdata.dipfit.rv;
end

%% stats
% pvals = std_stat({lsm_ersp_selfBase_pChange_mean; zeros(size(lsm_ersp_selfBase_pChange_mean))}, ...
%     'mode', 'fieldtrip', 'fieldtripmethod', 'montecarlo', 'condstats', 'on', 'fieldtripmcorrect', 'cluster');

[stats,df,pvals]= statcondfieldtrip({lsm_ersp_selfBase_pChange_mean; zeros(size(lsm_ersp_selfBase_pChange_mean))},...
    'paired','on',...
    'method','permutation',...
    'naccu', 2000,...
    'alpha', 0.05,...
    'mcorrect','cluster');

% [stats,df,pvals]= statcondfieldtrip({lsm_ersp_allBase_pChange_mean; zeros(size(lsm_ersp_allBase_pChange_mean))},...
%     'paired','on',...
%     'method','permutation',...
%     'naccu', 1000,...
%     'alpha', 0.05,...
%     'mcorrect','cluster');

[lsm_t_values, lsm_mask] = apply_significance_mask(pvals, stats);


[stats,df,pvals]= statcondfieldtrip({rsm_ersp_selfBase_pChange_mean; zeros(size(rsm_ersp_selfBase_pChange_mean))},...
    'paired','on',...
    'method','permutation',...
    'naccu', 2000,...
    'alpha', 0.05,...
    'mcorrect','cluster');

% [stats,df,pvals]= statcondfieldtrip({rsm_ersp_allBase_pChange_mean; zeros(size(rsm_ersp_allBase_pChange_mean))},...
%     'paired','on',...
%     'method','permutation',...
%     'naccu', 1000,...
%     'alpha', 0.05,...
%     'mcorrect','cluster');

[rsm_t_values, rsm_mask] = apply_significance_mask(pvals, stats);

saveName = fullfile(mainDir,'data',"study_p-masks.mat");
save(saveName,"rsm_mask","lsm_mask");

%% plot
% set parameters
colorRSM1 = {[0.90, 0.60, 0.00] * 0.9};
colorRSM2 = {[0.90, 0.60, 0.00]};
colorLSM1 = {[0.35, 0.70, 0.90] * 0.9};
colorLSM2 = {[0.35, 0.70, 0.90]};

% make plot
f1 = figure('units','normalized','outerposition',[0 0 1 1]);

colorList = [repmat(colorRSM2, length(rsm_dip), 1); repmat(colorLSM2, length(lsm_dip), 1)];
dipplot([rsm_dip, lsm_dip], ...
    'gui','off',...
    'verbose','off', ...
    'spheres','on', ...
    'dipolelength',0,...
    'view', [45 45], ...
    'color',colorList,...
    'dipolesize', [25*ones(length(rsm_dip),1);25*ones(length(lsm_dip),1)]);

% colorList = [colorRSM1; repmat(colorRSM2, length(rsm_dip), 1); colorLSM1; repmat(colorLSM2, length(lsm_dip), 1)];
% dipplot([computecentroid(rsm_dip), rsm_dip, computecentroid(lsm_dip), lsm_dip], ...
%     'gui','off',...
%     'verbose','off', ...
%     'spheres','on', ...
%     'dipolelength',0,...
%     'view', [45 45], ...
%     'color',colorList,...
%     'dipolesize', [55;25*ones(length(rsm_dip),1);55;25*ones(length(lsm_dip),1)]);


set(gcf,'Color',[1 1 1])

%%
f2 = figure('units','normalized','outerposition',[0 0 1 .50]);
tl = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');

maxVal = [];

ax2(1) = nexttile(1);
%y_mean = squeeze(mean(lsm_ersp_allBase_pChange_mean,3));
y_mean = squeeze(mean(lsm_ersp_selfBase_pChange_mean,3));

p1 = contourf(tfdata(1).times./1000,tfdata(1).freqs,y_mean,40,'linecolor','none');hold on
contour(tfdata(1).times./1000,tfdata(1).freqs,lsm_mask,1,'-k','LineWidth',1);
set(gca, 'ydir', 'normal'); ylabel('Frequency (Hz)','FontSize',18);
xlim([0 tfdata(1).warpVals(end)./1000]); xticks(tfdata(1).warpVals./1000);
xline(tfdata(1).warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
xline(tfdata(1).warpVals([2,4])./1000,'LineWidth',2,'alpha',1);
ylim([4 50]); yticks(tfdata(1).freqs(find(ismember(tfdata(1).freqs, [4 8 13 30 50]))))
xticklabels({'RHS','LTO','LHS','RTO','RHS'});xtickangle(30);
xlabel('Gait Cycle')
yline(13,'--k',LineWidth=2)
yline(30,'--k',LineWidth=2)
maxVal = [maxVal max(max(abs(y_mean)))];
%title('Left Sensorimotor Cluster','FontSize',21,'Interpreter','none');

ax2(2) = nexttile(2);
%y_mean = squeeze(mean(rsm_ersp_allBase_pChange_mean,3));
y_mean = squeeze(mean(rsm_ersp_selfBase_pChange_mean,3));
p1 = contourf(tfdata(1).times./1000,tfdata(1).freqs,y_mean,40,'linecolor','none');hold on
contour(tfdata(1).times./1000,tfdata(1).freqs,rsm_mask,1,'-k','LineWidth',1);
set(gca, 'ydir', 'normal');
xlim([0 tfdata(1).warpVals(end)./1000]); xticks(tfdata(1).warpVals./1000);
xline(tfdata(1).warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
xline(tfdata(1).warpVals([2,4])./1000,'LineWidth',2,'alpha',1);
ylim([4 50]); yticks(tfdata(1).freqs(find(ismember(tfdata(1).freqs, [4 8 13 30 50]))))
xticklabels({'RHS','LTO','LHS','RTO','RHS'});xtickangle(30);
xlabel('Gait Cycle')

yline(13,'--k',LineWidth=2)
yline(30,'--k',LineWidth=2)
maxVal = [maxVal max(max(abs(y_mean)))];
%title('Right Sensorimotor Cluster','FontSize',21,'Interpreter','none');

linkaxes(ax2)
ax2(1).XAxis.FontSize = 18;
ax2(2).XAxis.FontSize = 18;
ax2(1).YAxis.FontSize = 18;
ax2(2).YAxis.FontSize = 18;
maxVal = floor(max(maxVal));
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = '% Change';
cb.FontSize = 18;
set([ax2(1) ax2(2)], 'Colormap', redblue(100) , 'CLim', [-maxVal maxVal]);
set(gcf,'Color',[1 1 1])

% add annotation
% addSwingAnnotation(ax2(1), tfdata(1).warpVals(2)/1000, tfdata(1).warpVals(3)/1000, 'Ipsilateral Swing');
% addSwingAnnotation(ax2(1), tfdata(1).warpVals(4)/1000, tfdata(1).warpVals(5)/1000, 'Contralateral Swing');
% addSwingAnnotation(ax2(2), tfdata(1).warpVals(2)/1000, tfdata(1).warpVals(3)/1000, 'Contralateral Swing');
% addSwingAnnotation(ax2(2), tfdata(1).warpVals(4)/1000, tfdata(1).warpVals(5)/1000, 'Ipsilateral Swing');

%%% save
%sgtitle('ERSP Across the Gait Cycle','Fontsize',23)
filename = 'ERSP_averaged_ALL_percent.png';
saveas(gcf,fullfile(figFolder,filename));
filename = 'ERSP_averaged_ALL_percent.fig';
saveas(gcf,fullfile(figFolder,filename));
%close all

%% extract

% check pos and neg clusters for LSM
lsm_t_pos = lsm_t_values;
lsm_t_pos(lsm_t_values <= 0) = nan;
lsm_t_neg = lsm_t_values;
lsm_t_neg(lsm_t_values >= 0) = nan;,
lsm_sum_pos = sum(lsm_t_pos, 'all', 'omitnan');
lsm_sum_neg = sum(lsm_t_neg, 'all', 'omitnan');

figure
t2 = tiledlayout(2,2,'TileSpacing','Compact');

nexttile;
contourf(tfdata(1).times./1000, tfdata(1).freqs, lsm_t_pos, 40, 'linecolor', 'none'); 
colormap(gca, 'hot');
title(sprintf('LSM Positive Clusters (\\SigmaT = %.2f)', lsm_sum_pos));
nexttile;
contourf(tfdata(1).times./1000, tfdata(1).freqs, lsm_t_neg, 40, 'linecolor', 'none'); 
colormap(gca, 'abyss');
title(sprintf('LSM Negative Clusters (\\SigmaT = %.2f)', lsm_sum_neg));

% check pos and neg clusters for LSM
rsm_t_pos = rsm_t_values;
rsm_t_pos(rsm_t_values <= 0) = nan;
rsm_t_neg = rsm_t_values;
rsm_t_neg(rsm_t_values >= 0) = nan;,
rsm_sum_pos = sum(rsm_t_pos, 'all', 'omitnan');
rsm_sum_neg = sum(rsm_t_neg, 'all', 'omitnan');

nexttile;
contourf(tfdata(1).times./1000, tfdata(1).freqs, rsm_t_pos, 40, 'linecolor', 'none'); 
colormap(gca, 'hot');
title(sprintf('RSM Positive Clusters (\\SigmaT = %.2f)', rsm_sum_pos));
nexttile;
contourf(tfdata(1).times./1000, tfdata(1).freqs, rsm_t_neg, 40, 'linecolor', 'none'); 
colormap(gca, 'abyss');
title(sprintf('RSM Negative Clusters (\\SigmaT = %.2f)', rsm_sum_neg));



%% functions
function addSwingAnnotation(ax, t1, t2, label)
    axPos = ax.Position;
    xlims = ax.XLim;
    x1 = (t1 - xlims(1)) / diff(xlims);
    x2 = (t2 - xlims(1)) / diff(xlims);
    left = axPos(1) + axPos(3) * x1;
    width = axPos(3) * (x2 - x1);
    if contains(label, 'contra', 'IgnoreCase', true)
        annotation('textbox', [left, axPos(2) + axPos(4)-0.005, width, 0.06], ...
        'String', label, 'HorizontalAlignment', 'center', ...
        'EdgeColor', 'black','BackgroundColor', 'black',...
        'FontSize', 15, 'Color','white',...
        'VerticalAlignment', 'middle');
    else
         annotation('textbox', [left, axPos(2) + axPos(4)-0.005, width, 0.06], ...
        'String', label, 'HorizontalAlignment', 'center', ...
        'EdgeColor', 'white','BackgroundColor', 'white',...
        'FontSize', 15, 'Color','black',...
        'VerticalAlignment', 'middle');
    end
end
