function ERSP = removeOutliersERSP(ERSP,figFolder)
%% ERSP outlier removal

%parameters for PCA method
PCA_method_on = false;
explained_thres = .95;
PCA_thres = 3;

%parameters for Burst method
Burst_thres1 = 3;
threshold_ratio = 0.1;
Burst_thres2 = 5; % this applies to time and freq outliers

tf = ERSP.data.*conj(ERSP.data);
mean_tf = mean(tf,4);
baseA = squeeze(mean(mean_tf,3));
% log_mean_tf = 10*log10(mean_tf./baseA);
% log_tf      = 10*log10(tf./baseA);
p_tf = (tf - baseA)./baseA * 100;
new_tf = {};

for di = 1:size(tf,1)

    % di_log_tf = squeeze(log_tf(di,:,:,:));
    di_tf = squeeze(p_tf(di,:,:,:));

    % 1. Burst-based method
    tff = squeeze(sum(di_tf,1));
    med_tff = median(tff, 2);
    std_tff = std(tff,0,2);
    threshold1_high = med_tff + (Burst_thres1 * std_tff);
    threshold1_low  = med_tff - (Burst_thres1 * std_tff);
    threshold2_high = med_tff + (Burst_thres2 * std_tff);
    threshold2_low  = med_tff - (Burst_thres2 * std_tff);

    keep = true(size(tff, 2), 1);
    for r = 1:size(tff, 2) % Iterate over trials
        trial = tff(:,r);
        
        %cond 1: if 5% of data is over 3std
        exceed1 = trial > threshold1_high | trial < threshold1_low;
        condition1 = sum(exceed1) > (threshold_ratio * size(tff, 1));

        %cond 2: any point exceeds 5×std
        exceed2 = trial > threshold2_high | trial < threshold2_low;
        condition2 = any(exceed2);

        % Remove if either condition is true
        if condition1 || condition2
            keep(r) = false;
        end
    end
    rm1 = find(~keep); % track removes
    keep_idx1 = find(keep);

    % 2. Frequency deviation method
    tfp = squeeze(mean(di_tf,2));
    med_tfp = median(tfp, 2);
    std_tfp = std(tfp,0,2);
    threshold3_high = med_tfp + (Burst_thres2 * std_tfp);
    threshold3_low  = med_tfp - (Burst_thres2 * std_tfp);

    keep2 = true(size(tfp, 2), 1);
    for r = 1:size(tfp, 2)
        trial = tfp(:,r);
        
        %cond 1: if 5% of data is over 3std
        exceed1 = trial > threshold3_high | trial < threshold3_low;
        condition1 = sum(exceed1) > 0;

        % Remove if either condition is true
        if condition1
            keep2(r) = false;
        end
    end
    rm2 = find(~keep2); % track removes
    keep_idx2 = find(keep2);
    
    % 3. PCA method
    if PCA_method_on
        di_tf_3 = di_tf(:,:,keep);
        num_trials = size(di_tf_3, 3);
        data_matrix = reshape(di_tf_3, [], num_trials)'; % rows by features (aka each data point in TF)
        [coeff, score, latent] = pca(data_matrix);
        cumulative_variance = cumsum(latent) / sum(latent);
        num_components = find(cumulative_variance >= explained_thres, 1);
        reconstructed = score(:, 1:num_components) * coeff(:, 1:num_components)';
        reconstruction_error = sum((data_matrix - reconstructed).^2, 2);
        outlier_threshold = median(reconstruction_error) + PCA_thres * mad(reconstruction_error, 1);
        outlier_trials = reconstruction_error > outlier_threshold;
        rm3 = keep_idx(outlier_trials); % relative to the orginal indexes of di_tf
    else
        rm3 = [];
    end

    % combined removed index
    combined_rm = unique([rm1; rm2]);
    keepFinal = true(size(tff, 2), 1);
    keepFinal(combined_rm) = false;
    ERSP.keep{di} = keepFinal;
    ERSP.removedTrialNum{di} = length(combined_rm);
    new_tf(di) = {squeeze(tf(di,:,:,keepFinal))};

    %% plot outcomes
    f1 = figure('units','normalized','outerposition',[.1 .1 .8 .8]); hold on
    t1 = tiledlayout(2,1,'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile;
    if ~isempty(tff(:,~keep))
        p2 = plot(ERSP.times/1000,tff(:,~keep),'LineWidth',1); hold on
        set(p2, 'Color', [1 0 0 0.8]);
    end
    p1 = plot(ERSP.times/1000,tff(:,keep)); 
    set(p1, 'Color', [0 .5 1]);
    plot(ERSP.times/1000,median(tff,2),'k','LineWidth',5); hold on
    plot(ERSP.times/1000,threshold1_low,'--k','LineWidth',1.5);
    plot(ERSP.times/1000,threshold1_high,'--k','LineWidth',1.5);
    plot(ERSP.times/1000,threshold2_low,'--r','LineWidth',1.5);
    plot(ERSP.times/1000,threshold2_high,'--r','LineWidth',1.5);
    a_pos = ERSP.times(end)/1000;
    text(a_pos, threshold1_low(end), ' Lower Threshold 1', 'Color', 'k', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 10);
    text(a_pos, threshold1_high(end), ' Upper Threshold 1', 'Color', 'k', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
    text(a_pos, threshold2_low(end), ' Lower Threshold 2', 'Color', 'k', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 10);
    text(a_pos, threshold2_high(end), ' Upper Threshold 2', 'Color', 'k', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);
    ylabel('percent change from baseline')
    xlabel('gait cycle')
    title(['Temporal-based Removal. N = ' num2str(length(rm1)) '/' num2str(length(keep))])
    
    ax2 = nexttile;
    if ~isempty(rm2)
        p2 = plot(ERSP.freqs,tfp(:,~keep2),'LineWidth',1); hold on
        set(p2, 'Color', [1 0 0 0.8]);
    end
    p1 = plot(ERSP.freqs,tfp(:,keep2)); hold on
    set(p1, 'Color', [0 .5 1]);
    plot(ERSP.freqs,median(tfp, 2),'k','LineWidth',5); hold on
    plot(ERSP.freqs,threshold3_low,'--k','LineWidth',1.5);
    plot(ERSP.freqs,threshold3_high,'--k','LineWidth',1.5);
    ylabel('average Power')
    xlabel('frequency')
    title(['Frequency-based Removal. N = ' num2str(length(rm2)) '/' num2str(length(keep2))])

    if ~isempty(rm3)
        ax3 = nexttile;
        x = 1:length(reconstruction_error);
        p3 = bar(x, reconstruction_error, 'BarWidth', 0.8);hold on;
        p3.FaceColor = 'flat';
        p3.CData(outlier_trials, :) = repmat([1 0 0], length(rm3), 1);
        yline(outlier_threshold, '--r', 'LineWidth', 2);
        ylabel('reconstruction error')
        title(['PCA Method Outlier Removal N = ' num2str(length(rm3)) '/' num2str(length(keep))])
    end

    set(gcf,'Color',[1 1 1])
    saveName = ['trialsRemoved-' ERSP.subject '_dipfit-' num2str(di) '.png'];
    saveas(gcf,fullfile(figFolder,saveName));
    close all

    %% plot ERSP results
    f2 = figure('units','normalized','outerposition',[.1 .1 .8 .8]);
    t0=tiledlayout(2,2);
    ymax = [];
    sgtitle([ERSP.subject '  Dipole-' num2str(di) '  Area-' ERSP.dipfit.model(di).areadk  '  RV-' num2str(ERSP.dipfit.model(di).rv)],...
         'Interpreter', 'none')

    t1 = tiledlayout(t0,2,1,'TileSpacing', 'compact', 'Padding', 'compact');
    t1.Layout.Tile=1;
    t1.Layout.TileSpan=[2 1];
    
    ax1 = nexttile(t1);
    y_mean = squeeze(mean(di_tf(:,:,rm1),3));
    contourf(ERSP.times/1000, ERSP.freqs, y_mean,40,'linecolor','none');
    ymax = [ymax max(max(y_mean))];
    xlim([0 ERSP.warpVals(end)./1000]); xticks(ERSP.warpVals./1000);
    xline(ERSP.warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
    xline(ERSP.warpVals([2,4])./1000,'LineWidth',2,'alpha',1);
    xticklabels({'RHS','LTO','LHS','RTO','RHS'}); xtickangle(45)
    ylim([2 55]); yticks(ERSP.freqs(find(ismember(ERSP.freqs, [8 13 30]))))
    title(['Temporal-based Removal. N = ' num2str(length(rm1)) ' / ' num2str(length(keep))])

    ax2 = nexttile(t1);
    if ~isempty(rm2)
        y_mean = squeeze(mean(di_tf(:,:,rm2),3));
        contourf(ERSP.times/1000, ERSP.freqs, y_mean,40,'linecolor','none');
        ymax = [ymax max(max(y_mean))];
        xlim([0 ERSP.warpVals(end)./1000]); xticks(ERSP.warpVals./1000);
        xline(ERSP.warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
        xline(ERSP.warpVals([2,4])./1000,'LineWidth',2,'alpha',1);
        xticklabels({'RHS','LTO','LHS','RTO','RHS'}); xtickangle(45)
        ylim([2 55]); yticks(ERSP.freqs(find(ismember(ERSP.freqs, [8 13 30]))))
        title(['Frequency-Based Removal N = ' num2str(length(rm2)) ' / ' num2str(length(keep))])
    end

    % ax3 = nexttile(t1);
    % y_mean = squeeze(mean(di_tf(:,:,keepFinal),3));
    % contourf(ERSP.times/1000, ERSP.freqs, y_mean,40,'linecolor','none');
    % xlim([0 ERSP.warpVals(end)./1000]); xticks(ERSP.warpVals./1000);
    % xline(ERSP.warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
    % xline(ERSP.warpVals([2,4])./1000,'LineWidth',2,'alpha',1);
    % xticklabels({'RHS','LTO','LHS','RTO','RHS'}); xtickangle(45)
    % ylim([2 55]); yticks(ERSP.freqs(find(ismember(ERSP.freqs, [8 13 30]))))
    % title(['Remaining Cycles. N = ' num2str(length(find(keepFinal == true))) ' / ' num2str(length(keep))])
    
    ymax = max(squeeze([ymax max(max(y_mean))]));
    set([ax1 ax2], 'Colormap', redblue(100), 'CLim', [-ymax ymax], 'YDir', 'normal');
    linkaxes([ax1 ax2])
    cb = colorbar(ax1);
    cb.Layout.Tile = 'east';
    cb.Label.String = '%';

    t2 = tiledlayout(t0,2,1,'TileSpacing', 'compact', 'Padding', 'compact');
    t2.Layout.Tile=2;
    t2.Layout.TileSpan=[2 1];
    ymax = [];
    
    ax4 = nexttile(t2);
    y_mean = squeeze(mean(di_tf,3));
    contourf(ERSP.times/1000, ERSP.freqs, y_mean,40,'linecolor','none');
    ymax = max(squeeze([ymax max(max(y_mean))]));
    xlim([0 ERSP.warpVals(end)./1000]); xticks(ERSP.warpVals./1000);
    xline(ERSP.warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
    xline(ERSP.warpVals([2,4])./1000,'LineWidth',2,'alpha',1);
    xticklabels({'RHS','LTO','LHS','RTO','RHS'}); xtickangle(45)
    ylim([2 55]); yticks(ERSP.freqs(find(ismember(ERSP.freqs, [8 13 30]))))
    title(['Original ERSP. N = ' num2str(length(keep)) ' / ' num2str(length(keep))])
    ymax = [ymax max(max(y_mean))];

    ax5 = nexttile(t2);
    %%% remake tf using new baseline
    new_tf2plot = squeeze(tf(di,:,:,keepFinal));
    mean_new_tf2plot = mean(new_tf2plot,3);
    new_baseA2plot = squeeze(mean(mean_new_tf2plot,2));
    di_new_tf_2plot = (new_tf2plot - new_baseA2plot)./new_baseA2plot *100;
    %%%
    y_mean = squeeze(mean(di_new_tf_2plot,3));
    contourf(ERSP.times/1000, ERSP.freqs, y_mean,40,'linecolor','none');
    ymax = max(squeeze([ymax max(max(y_mean))]));
    xlim([0 ERSP.warpVals(end)./1000]); xticks(ERSP.warpVals./1000);
    xline(ERSP.warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
    xline(ERSP.warpVals([2,4])./1000,'LineWidth',2,'alpha',1);
    xticklabels({'RHS','LTO','LHS','RTO','RHS'}); xtickangle(45)
    ylim([2 55]); yticks(ERSP.freqs(find(ismember(ERSP.freqs, [8 13 30]))))
    title(['Cleaned ERSP. N = ' num2str(length(find(keepFinal == true))) ' / ' num2str(length(keep))])
    ymax = max([ymax max(max(y_mean))]);
    
    set([ax4 ax5], 'Colormap', redblue(100), 'CLim', [-ymax ymax], 'YDir', 'normal');
    linkaxes([ax4 ax5])
    cb = colorbar(ax4);
    cb.Layout.Tile = 'east';
    cb.Label.String = '%';
    
    set(gcf,'Color',[1 1 1]);
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 17);
    saveName = ['erspRemoved-' ERSP.subject '_dipole-' num2str(di) '.fig'];
    saveas(gcf,fullfile(figFolder,saveName));
    %close all

    % visualize
    %figure
    % for i = 1:length(rm1)
    %     nexttile
    %     imagesc(ERSP.times/1000, ERSP.freqs, di_tf(:,:,rm1(i)));
    %     set(gca, 'YDir', 'normal'); % Ensure correct frequency ordering
    %     colormap jet;
    %     colorbar;
    % end
    % 
    % for i = 1:length(rm2)
    %     nexttile
    %     imagesc(ERSP.times/1000, ERSP.freqs, di_tf(:,:,rm2(i)));
    %     set(gca, 'YDir', 'normal'); % Ensure correct frequency ordering
    %     colormap jet;
    %     colorbar;
    % end
    % for r = 1:length(combined_rm)
    %     nexttile
    %     ymax = max(max(di_tf(:,:,combined_rm(r))));
    %     contourf(ERSP.times/1000, ERSP.freqs, di_tf(:,:,combined_rm(r)),40,'linecolor','none');
    %     set(gca, 'YDir', 'normal'); % Ensure correct frequency ordering
    %     set(gca, 'Colormap', jet, 'CLim', [-ymax ymax]);
    %     colorbar;
    % end
end

ERSP.data_cleaned = new_tf;
end