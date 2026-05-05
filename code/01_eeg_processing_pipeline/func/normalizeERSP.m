function tfdata = normalizeERSP(ERSP)
tfdata = [];

% update gait table if ERSP has been cleaned
if isfield(ERSP,'data_cleaned')
    ERSP.gait_table = ERSP.gait_table(ERSP.keep{idx_select},:);
end

% get condition indexes
idx_nocue = find(strcmp(ERSP.gait_table.Condition,'walk-run1'));
if isempty(idx_nocue)
    idx_nocue = find(strcmp(ERSP.gait_table.Condition,'walk-run2'));
end
idx_audi = find(strcmp(ERSP.gait_table.Condition,'walkAuditory-run1'));
if isempty(idx_audi)
    idx_audi = find(strcmp(ERSP.gait_table.Condition,'walkAuditory-run2'));
end
idx_visu = find(strcmp(ERSP.gait_table.Condition,'walkVisual-run1'));
if isempty(idx_visu)
    idx_visu = find(strcmp(ERSP.gait_table.Condition,'walkVisual-run2'));
end

tf = ERSP.data.*conj(ERSP.data);
tf_mean = squeeze(mean(tf,4,"omitnan"));

for idx_select = 1:size(ERSP.data,1)

tmp_data_nocue = []; tmp_data_audi = []; tmp_data_visu = [];

tf_select = squeeze(tf_mean(idx_select,:,:));
baseA = squeeze(ERSP.baselineAllCond(idx_select,:))';

mean_db_data = 10*log10(bsxfun(@rdivide, tf_select, baseA));
mean_p_data = (tf_select - baseA)./baseA *100;

tmp_data_nocue = squeeze(tf(idx_select,:,:,idx_nocue));
tmp_data_audi = squeeze(tf(idx_select,:,:,idx_audi));
tmp_data_visu = squeeze(tf(idx_select,:,:,idx_visu));

%%%%% ersp by condition %%%%%

% no cue %%%
tf_nocue_mean = squeeze(mean(tmp_data_nocue,3,"omitnan"));
tf_nocue_base = mean(tf_nocue_mean,2,"omitnan");

ersp_allBase_pChange{1} = (tmp_data_nocue - baseA)./baseA * 100;
ersp_allBase_db_mean{1}      = 10*log10(bsxfun(@rdivide, tf_nocue_mean, baseA));

ersp_selfBase_pChange{1}    = (tmp_data_nocue - tf_nocue_base)./tf_nocue_base * 100;
ersp_selfBase_db_mean{1}         = 10*log10(bsxfun(@rdivide, tf_nocue_mean, tf_nocue_base));

gait_table{1}   = ERSP.gait_table(idx_nocue,:);

% audi %%%
tf_audi_mean = squeeze(mean(tmp_data_audi,3,"omitnan"));
tf_audi_base = mean(tf_audi_mean,2,"omitnan");

ersp_allBase_pChange{2} = (tmp_data_audi - baseA)./baseA * 100;
ersp_allBase_db_mean{2}      = 10*log10(bsxfun(@rdivide, tf_audi_mean, baseA));

ersp_selfBase_pChange{2} = (tmp_data_audi - tf_audi_base)./tf_audi_base * 100;
ersp_selfBase_db_mean{2}      = 10*log10(bsxfun(@rdivide, tf_audi_mean, tf_audi_base));

ersp_mean_diff_p_allBase{1}                 = mean(ersp_allBase_pChange{2},3) - mean(ersp_allBase_pChange{1},3);
ersp_mean_diff_db_allBase{1}                = mean(ersp_allBase_db_mean{2},3) - mean(ersp_allBase_db_mean{1},3);

ersp_mean_diff_p_selfBase{1}                 = mean(ersp_selfBase_pChange{2},3) - mean(ersp_selfBase_pChange{1},3);
ersp_mean_diff_db_selfBase{1}                = mean(ersp_selfBase_db_mean{2},3) - mean(ersp_selfBase_db_mean{1},3);

gait_table{2}   = ERSP.gait_table(idx_audi,:);

% visual %%%
tf_visu_mean = squeeze(mean(tmp_data_visu,3,"omitnan"));
tf_visu_base = mean(tf_visu_mean,2,"omitnan");

ersp_allBase_pChange{3} = (tmp_data_visu - baseA)./baseA * 100;
ersp_allBase_db_mean{3}      = 10*log10(bsxfun(@rdivide, tf_visu_mean, baseA));

ersp_selfBase_pChange{3}    = (tmp_data_visu - tf_visu_base)./tf_visu_base * 100;
ersp_selfBase_db_mean{3}         = 10*log10(bsxfun(@rdivide, tf_visu_mean, tf_visu_base));

ersp_mean_diff_p_allBase{2}                 = mean(ersp_allBase_pChange{3},3) - mean(ersp_allBase_pChange{1},3);
ersp_mean_diff_db_allBase{2}                = mean(ersp_allBase_db_mean{3},3) - mean(ersp_allBase_db_mean{1},3);

ersp_mean_diff_p_selfBase{2}                 = mean(ersp_selfBase_pChange{3},3) - mean(ersp_selfBase_pChange{1},3);
ersp_mean_diff_db_selfBase{2}                = mean(ersp_selfBase_db_mean{3},3) - mean(ersp_selfBase_db_mean{1},3);

gait_table{3}   = ERSP.gait_table(idx_visu,:);

% save to struct
tfdata(idx_select).subject                  = ERSP.subject;
tfdata(idx_select).group                    = ERSP.group;
tfdata(idx_select).HS_idx                   = ERSP.HS_idx;

tfdata(idx_select).ersp_selfBase_pChange        = ersp_selfBase_pChange;
tfdata(idx_select).ersp_selfBase_pChange_mean   = cellfun(@(x) mean(x, 3), ersp_selfBase_pChange, 'UniformOutput', false);

tfdata(idx_select).ersp_allBase_pChange         = ersp_allBase_pChange;
tfdata(idx_select).ersp_allBase_pChange_mean    = cellfun(@(x) mean(x, 3), ersp_allBase_pChange, 'UniformOutput', false);

tfdata(idx_select).ersp_selfBase_db_mean        = ersp_selfBase_db_mean;
tfdata(idx_select).ersp_allBase_db_mean         = ersp_allBase_db_mean;

tfdata(idx_select).ersp_db_allCond_mean         = mean_db_data; % average across all cycles regardless of condition
tfdata(idx_select).ersp_p_allCond_mean          = mean_p_data;

tfdata(idx_select).ersp_diff_pChange_allBase    = ersp_mean_diff_p_allBase;
tfdata(idx_select).ersp_diff_db_allBase         = ersp_mean_diff_db_allBase;

tfdata(idx_select).ersp_diff_pChange_selfBase   = ersp_mean_diff_p_selfBase;
tfdata(idx_select).ersp_diff_db_selfBase        = ersp_mean_diff_db_selfBase;

tfdata(idx_select).labels                   = {'nocue' 'auditory' 'visual'};
tfdata(idx_select).diff_labels              = {'auditory - nocue' 'visual - nocue'};
tfdata(idx_select).cond_index               = {idx_nocue,idx_audi,idx_visu};
tfdata(idx_select).gait_table               = gait_table;
tfdata(idx_select).times                    = ERSP.times;
tfdata(idx_select).freqs                    = ERSP.freqs;
tfdata(idx_select).warpVals                 = ERSP.warpVals;
tfdata(idx_select).dipfit           = ERSP.dipfit.model(idx_select);
tfdata(idx_select).chanlocs         = ERSP.chanlocs;
tfdata(idx_select).icawinv          = ERSP.icawinv(:,idx_select);

end