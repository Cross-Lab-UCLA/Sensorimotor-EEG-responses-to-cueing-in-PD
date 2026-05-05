function plotERSP(ERSP,ersp_all,f_bound,min_idx_list,idx,figFolder,roi_name)

ylim_vals = [4 50];

if ~isfield(ERSP,'tfdata')
    disp('Time frequency data does not exist in current ERSP set!')
    keyboard
end
%%
figure('units','normalized','outerposition',[0 0 1 1])
set(gcf,'Color',[1 1 1])
maxVal = [];

if ndims(ersp_all) == 3
    ersp_all = permute(ersp_all,[2 3 1]);
end

tl = tiledlayout(size(ersp_all,3)*2,12, 'TileSpacing', 'tight', 'Padding', 'compact');

for p = 1:size(ersp_all,3)

    ax1(p) = nexttile(tl);
    ax1(p).Layout.TileSpan = [2 2];
    dip = ERSP.dipfit.model(min_idx_list(p));
    p1 = dipplot(dip, 'coordformat', 'MNI', ...
        'gui', 'off', 'dipolesize', 30, 'view', [45 35],...
        'projlines', 'on','verbose','off'); hold on
    axis(ax1(p), 'image', 'tight');

    ax2(p) = nexttile(tl);
    ax2(p).Layout.TileSpan = [2 2];
    axes(ax2(p));
    dip = ERSP.dipfit.model(min_idx_list(p));
    p2 = topoplot(ERSP.icawinv(:,min_idx_list(p)),ERSP.chanlocs);
    axis(ax2(p), 'image', 'tight'); 

    ax3(p) = nexttile(tl);
    ax3(p).Layout.TileSpan = [2 4];
    data = squeeze(ersp_all(:,:,p));
    contourf(ERSP.times/1000,ERSP.freqs,data,40,'linecolor','none');hold on
    set(gca, 'ydir', 'normal');
    ylabel('Frequency (Hz)','FontSize',13);
    xlim([0 ERSP.warpVals(end)./1000]); xticks(ERSP.warpVals./1000);
    xline(ERSP.warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
    xline(ERSP.warpVals([2,4])./1000,'LineWidth',1,'alpha',1);
    yline(ERSP.freqs(f_bound),'--r');
    ylim(ylim_vals); yticks(ERSP.freqs(find(ismember(ERSP.freqs, [8 13 30]))))
    xticklabels({'RHS','LTO','LHS','RTO','RHS'});
    xtickangle(45)
    maxVal = max([maxVal max(max(abs(data)))]);
    title({dip.areadk [' rv-' num2str(dip.rv)]},'Interpreter','none');
    cb = colorbar(ax2(p), 'Location', 'eastoutside');
    cb.Label.String = '%';
    colormap(redblue(100));

    ax4(p) = nexttile(tl);
    ax4(p).Layout.TileSpan = [2 4];
    y2_mean = squeeze(mean(ersp_all(f_bound(1):f_bound(2),:,p),1));
    plot(ERSP.times/1000,y2_mean,'LineWidth',3.5,'Color','k');hold on
    yline(0,'--r');
    set(gca, 'ydir', 'normal');
    xlim([0 ERSP.warpVals(end)./1000]); xticks(ERSP.warpVals./1000);
    xline(ERSP.warpVals([1,3,5])./1000,'LineWidth',2,'alpha',1);
    xline(ERSP.warpVals([2,4])./1000,'LineWidth',1,'alpha',1);
    xtickangle(45)
    xticklabels({'RHS','LTO','LHS','RTO','RHS'});
    ylabel('Beta Power (%)','FontSize',13);


    if p == idx
        ax3(p).XColor = 'r';
        ax3(p).YColor = 'r';
        ax3(p).LineWidth = 3;
        ax4(p).XColor = 'r';
        ax4(p).YColor = 'r';
        ax4(p).LineWidth = 3;
    end

end


set(findall(ax3, '-property', 'FontSize'), 'FontSize', 20)
set(findall(ax4, '-property', 'FontSize'), 'FontSize', 20)
linkaxes(ax3)
linkaxes(ax4)
set(ax3, 'CLim', [-maxVal maxVal]);
set(gcf,'Color',[1 1 1])
tname = [ERSP.subject '_' roi_name];
sgtitle(tname,'Interpreter', 'none');

%% Save

folderName = fullfile(figFolder,['ERSP_' roi_name]);
if ~exist(folderName, 'dir')
    mkdir(folderName)
end
saveas(gcf,fullfile(folderName,[tname '.png']))
saveas(gcf,fullfile(folderName,[tname '.fig']))
close all
end