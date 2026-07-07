function [data_group_map,data_topo_corr,data_icawinv_norm,data_subject] = getTopoCorr(data,nBoot,useOldMedian)
%% getTopoCorr
%
% input:
% data = data struct containing [icawinv, subject] for all trials or ersps
%        selected for analysis
% nBoot = number of iterations for picking representative subject icawinv.
%         If nBoot = [] or 1, no bootstrapping is performed and the group
%         map is a subject-balanced median (one representative IC per subject).
% useOldMedian = optional logical (default false). If true and not
%         bootstrapping, use the old method: coordinate-wise median across all
%         pooled ICs. 
%
% output:
% data_group_map    = group centrotype icawinv
% data_topo_corr    = corr r value of each icawinv to the group map
% data_icawinv_norm = RMS-normalized icawinv values
% data_subject      = subject label for each column of data_icawinv_norm /
%                      each entry of data_topo_corr (same order as input data)
%
% LM 071426
%%
if nargin < 3 || isempty(useOldMedian)
    useOldMedian = false;
end

data_icawinv = [data(:).icawinv];
data_subject = {data.subject};

rms_winv = rms(data_icawinv);
data_icawinv_norm = data_icawinv ./ rms_winv;
subj_names = unique(data_subject);

nChan = size(data_icawinv_norm, 1);

if ~isempty(nBoot) & nBoot > 1
    disp(['running bootstrap method...iter @ ' num2str(nBoot)])
    rs = RandStream('mt19937ar','Seed',5489);
    boot_maps = zeros(nChan, nBoot);
    for b = 1:nBoot
        sel = zeros(nChan, length(subj_names));
        for s = 1:length(subj_names)
            idx = find(strcmp(data_subject, subj_names{s}));
            pick = idx(randi(rs, numel(idx)));
            sel(:, s) = data_icawinv_norm(:, pick);
        end
        boot_maps(:, b) = mean(sel, 2);
    end
    data_group_map = median(boot_maps, 2);
elseif useOldMedian
    % old method: coordinate-wise median across all pooled ICs
    data_group_map = median(data_icawinv_norm, 2);
else
    % deterministic subject-balanced reference (each subject weighted equally)
    R_all = corr(data_icawinv_norm, 'rows','pairwise');

    subj_maps     = nan(nChan, numel(subj_names));
    subj_best_idx = nan(numel(subj_names),1);
    for s = 1:numel(subj_names)
        mask = strcmp(data_subject, subj_names{s});
        is   = find(mask);
        otherNames = subj_names(~strcmp(subj_names, subj_names{s}));
        cand_sim = nan(numel(is), numel(otherNames));
        for o = 1:numel(otherNames)
            oIdx = strcmp(data_subject, otherNames{o});
            cand_sim(:,o) = max(R_all(is, oIdx), [], 2, 'omitnan');   % that other subject's single best-matching candidate
        end
        xSubjSim = mean(cand_sim, 2, 'omitnan'); % one balanced vote per other subject
        [~, best] = max(xSubjSim);
        subj_best_idx(s) = is(best);
        subj_maps(:,s)   = data_icawinv_norm(:, is(best));
    end

    % the subject representative most correlated with all others
    R = R_all(subj_best_idx, subj_best_idx);
    R(1:size(R,1)+1:end) = 0;
    [~, ctrIdx] = max(sum(R, 2, 'omitnan'));
    data_group_map = subj_maps(:, ctrIdx);
end

data_topo_corr = corr(data_icawinv_norm, data_group_map, 'rows','pairwise')';
