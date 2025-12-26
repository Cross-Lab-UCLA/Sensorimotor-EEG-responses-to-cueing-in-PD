function [data_group_map data_topo_corr,data_icawinv_norm] = getTopoCorr(data,nBoot)
%%
%
% input:
% data = data strcut containing [icawinv, subject] for all trials or ersps
% selected for analysis
% 
% nBoot = number of iterations for picking repersentaive subject icawinv
%         if nBoot = [] or 1, no bootstrapping is performed. 
% output:
% data_group_map    = group median icawinv
% data_topo_corr    = corr r values for each icawinv
% data_icawinv_norm = normalized icawinv values
%
% LM 012825
%%
data_icawinv = [data(:).icawinv];

rms_winv = rms(data_icawinv);
data_icawinv_norm = data_icawinv ./ rms_winv;
subj_names = unique({data.subject});

[nChan, nIC] = size(data_icawinv_norm);

if ~isempty(nBoot) & nBoot > 1
    disp(['running bootstrap method...iter @ ' num2str(nBoot)])
    boot_maps = zeros(nChan, nBoot);
    for b = 1:nBoot
        sel = zeros(nChan, length(subj_names));
        for s = 1:length(subj_names)
            idx = find(strcmp({data.subject}, subj_names{s}));
            pick = idx(randi(numel(idx)));
            sel(:, s) = data_icawinv_norm(:, pick);
        end
        boot_maps(:, b) = mean(sel, 2);
    end
    data_group_map = median(boot_maps, 2);
    data_topo_corr = zeros(1, nIC);
    for i = 1:nIC
        data_topo_corr(i) = corr(data_icawinv_norm(:,i), data_group_map, 'rows','pairwise');
    end

else
    data_group_map = median(data_icawinv_norm, 2);
    for i = 1:nIC
        data_topo_corr(i) = corr(data_icawinv_norm(:,i), data_group_map, 'rows','pairwise');
    end
end