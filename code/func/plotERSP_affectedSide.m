function plotERSP_affectedSide(ERSP,ersp_all,f_bound,min_idx_list,idx,figFolder,roi_name,affectedSide)

if ~isfield(ERSP,'data')
    disp('Time frequency data does not exist in current ERSP set!')
    keyboard
end

figure('units','normalized','outerposition',[0 0 1 1])
set(gcf,'Color',[1 1 1])
maxVal = [];

if ndims(ersp_all) == 3
    ersp_all = permute(ersp_all,[2 3 1]);
end

tl = tiledlayout(size(ersp_all,3),8, 'TileSpacing', 'tight', 'Padding', 'compact');

for p = 1:size(ersp_all,3)

    ax1(p) = nexttile;
    dip = ERSP.dipfit.model(min_idx_list(p));
    p1 = dipplot(dip, 'coordformat', 'MNI', ...
        'gui', 'off', 'dipolesize', 30, 'view', [45 35],...
        'projlines', 'on','verbose','off'); hold on
    set(gcf,'Color','w');

    ax2(p) = nexttile([1 3]);
    data = squeeze(ersp_all(:,:,p));
    contourf(ERSP.times/1000,ERSP.freqs,data,40,'linecolor','none');hold on
    set(gca, 'ydir', 'normal'); ylabel('Frequencies (Hz)','FontSize',11);
    xlim([0 ERSP.warpVals(end)./1000]); xticks(ERSP.warpVals./1000);
    xline(ERSP.warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
    xline(ERSP.warpVals([2,4])./1000,'LineWidth',1,'alpha',1);
    yline(ERSP.freqs(f_bound),'--r');
    ylim([2 55]); yticks(ERSP.freqs(find(ismember(ERSP.freqs, [8 13 30]))))
    xticklabels({'MoreAffected HS','LessAffected TO','LessAffect HS','MoreAffected TO','MoreAffected HS'});
    xtickangle(45)
    maxVal = max([maxVal max(max(abs(data)))]);
    title({dip.areadk [' rv-' num2str(dip.rv)]},'Interpreter','none');
    cb = colorbar(ax2(p), 'Location', 'eastoutside');
    cb.Label.String = '%';
    colormap(redblue(100));

    ax3(p) = nexttile([1 3]);
    y2_mean = squeeze(mean(ersp_all(f_bound(1):f_bound(2),:,p),1));
    plot(ERSP.times/1000,y2_mean);hold on
    yline(0,'--r');
    set(gca, 'ydir', 'normal');
    xlim([0 ERSP.warpVals(end)./1000]); xticks(ERSP.warpVals./1000);
    xline(ERSP.warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
    xline(ERSP.warpVals([2,4])./1000,'LineWidth',1,'alpha',1);
    xtickangle(45)
    xticklabels({'MoreAffected HS','LessAffected TO','LessAffect HS','MoreAffected TO','MoreAffected HS'});
    title('mean(beta % change) across time','Interpreter','none');

    if p == idx
        ax2(p).XColor = 'r';
        ax2(p).YColor = 'r';
        ax2(p).LineWidth = 3;
        ax3(p).XColor = 'r';
        ax3(p).YColor = 'r';
        ax3(p).LineWidth = 3;
    end

end
set(findall(ax2, '-property', 'FontSize'), 'FontSize', 20)
set(findall(ax3, '-property', 'FontSize'), 'FontSize', 20)
linkaxes(ax3)
set(ax2, 'CLim', [-maxVal maxVal]);
linkaxes(ax1)
linkaxes(ax2)

tname = [ERSP.subject '_' affectedSide 'AffectedSide_' roi_name];
sgtitle(tname,'Interpreter', 'none');
folderName = fullfile(figFolder,['ERSP_' affectedSide]);
if ~exist(folderName, 'dir')
    mkdir(folderName)
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 24)

saveas(gcf,fullfile(folderName,[tname '.png']))
saveas(gcf,fullfile(folderName,[tname '.fig']))
close all
end