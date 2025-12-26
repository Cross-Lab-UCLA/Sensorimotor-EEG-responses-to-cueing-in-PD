function [group_pattern, data_ersp_vec,nFreq,nTime] = getERSPCorr(data,nBoot)
%%
%
% input:
% data = data strcut containing [icawinv, subject] for all trials or ersps
% selected for analysis
% nBoot = number of iterations for picking repersentaive subject icawinv
%
% output:
% data_group_map = group median icawinv
% data_topo_corr = corr r values for each icawinv
%%

data_ersp = cat(3, data(:).ersp_selfBase_pChange_mean);
data_ersp_vec = reshape(data_ersp, [], size(data_ersp,3));  % [nFeat x nICs]
subj_names = unique({data.subject});
nSubj = numel(subj_names);
nBoot = 2000;
nFreq = size(data_ersp,1);
nTime = size(data_ersp,2);

[nFeat, nIC] = size(data_ersp_vec);
boot_maps = zeros(nFeat, nBoot);

for b = 1:nBoot
    sel = zeros(nFeat, nSubj);
    for s = 1:nSubj
        idx = find(strcmp({data.subject}, subj_names{s}));
        pick = idx(randi(numel(idx)));      % pick ONE ERSP from this subject
        sel(:,s) = data_ersp_vec(:,pick);
    end
    boot_maps(:,b) = median(sel,2,'omitnan'); % subject-balanced mean ERSP
end

group_pattern = median(boot_maps, 2);
