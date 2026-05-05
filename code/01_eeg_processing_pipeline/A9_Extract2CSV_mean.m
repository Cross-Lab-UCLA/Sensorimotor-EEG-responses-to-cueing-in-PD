%% create csv from Targeted Clusters
%
% Run after A7_Target_Clusters. Takes RSM and LSM and extract relevant
% variables for statistical analysis.
%
% LM 072525
%%
clear all; clc;
eeglab; close all;
if ispc
    mainDir     = 'C:\Git\DoD-Gait';
else
    mainDir = '/Users/Leo/Git/DoD-Gait';
end
dataDir = fullfile(mainDir,'data');
saveDir = fullfile(mainDir,'results','processed_data');
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end
% prompt user for folder
%path = uigetdir(dataDir);
%label = extractAfter(path,'selected_');

path = fullfile(dataDir,'07_selected_withStanding_run1only_moreAffected');
clusterLabel = {'RSM' 'LSM'};

%% get mask
p_maskfile = fullfile(mainDir,'data',"study_p-masks.mat");
if exist(p_maskfile, 'file')
    load(p_maskfile)
end

%% run loop
for clu = 1:length(clusterLabel)
    Files = dir(fullfile(path,[clusterLabel{clu} '*']));
    if exist(p_maskfile, 'file') && strcmp(clusterLabel{clu}, 'RSM')
        p_mask = rsm_mask;
    elseif exist(p_maskfile, 'file') && strcmp(clusterLabel{clu}, 'LSM')
        p_mask = lsm_mask;
    else
        p_mask = [];
    end

    % organize data
    ersp_data = [];
    T = table();
    counter = 0;
    for s = 1:length(Files)

        disp(['Loading ' Files(s).name])
        load(fullfile(Files(s).folder,Files(s).name));
        disp('Loading complete.');

        curr_subj = tfdata.subject;
        disp(['Processing ' curr_subj]);
        times = tfdata.times;
        freqs = tfdata.freqs;
        warpVals = tfdata.warpVals;
        t_bound(1) = find(times == warpVals(1)); % RHS
        [~, t_bound(2)] = min(abs(times - warpVals(2))); % LTO
        [~, t_bound(3)] = min(abs(times - warpVals(3))); % LHS
        [~, t_bound(4)] = min(abs(times - warpVals(4))); % RTO
        [~, t_bound(5)] = min(abs(times - warpVals(5))); % RHS
        theta_f_bound = [find(freqs == 4) find(freqs == 12)];
        mu_f_bound = [find(freqs == 7) find(freqs == 12)];
        beta_f_bound = [find(freqs == 13) find(freqs == 30)];

        for c = 1:length(tfdata.labels)
            % get data
            data_allBase_db         = tfdata.ersp_allBase_db_mean{c};
            data_allBase_pChange    = tfdata.ersp_allBase_pChange_mean{c};
            data_selfBase_db        = tfdata.ersp_selfBase_db_mean{c};
            data_selfBase_pChange   = tfdata.ersp_selfBase_pChange_mean{c};

            if ~isempty(p_mask)
                p_mask(p_mask == 0) = NaN;
                data_allBase_pChange_masked = data_allBase_pChange.*p_mask;
                beta_allBase_mean_pChange_masked       = squeeze(mean(data_allBase_pChange_masked(beta_f_bound(1):beta_f_bound(2), :),'omitmissing'));
                beta_allBase_mean_pChange_masked_left_swing   = beta_allBase_mean_pChange_masked(t_bound(2):t_bound(3));
                beta_allBase_mean_pChange_masked_right_swing  = beta_allBase_mean_pChange_masked(t_bound(4):t_bound(5));

                data_selfBase_pChange_masked = data_selfBase_pChange.*p_mask;
                beta_selfBase_mean_pChange_masked       = squeeze(mean(data_selfBase_pChange_masked(beta_f_bound(1):beta_f_bound(2), :),'omitmissing'));
                beta_selfBase_mean_pChange_masked_left_swing   = beta_selfBase_mean_pChange_masked(t_bound(2):t_bound(3));
                beta_selfBase_mean_pChange_masked_right_swing  = beta_selfBase_mean_pChange_masked(t_bound(4):t_bound(5));
            end

            % Compute allBase means for dB data
            beta_allBase_mean_db                    = squeeze(mean(data_allBase_db(beta_f_bound(1):beta_f_bound(2), :)));
            beta_allBase_mean_db_left_swing         = beta_allBase_mean_db(t_bound(2):t_bound(3));
            beta_allBase_mean_db_leftMid2Stance     = beta_allBase_mean_db(floor((t_bound(2) + t_bound(3))/2):t_bound(4)-1);
            beta_allBase_mean_db_right_swing        = beta_allBase_mean_db(t_bound(4):t_bound(5));
            beta_allBase_mean_db_rightMid2Stance    = [ beta_allBase_mean_db(floor((t_bound(4) + t_bound(5))/2):t_bound(5)) ...
                beta_allBase_mean_db(t_bound(1):t_bound(2)-1)];

            mu_allBase_mean_db                = squeeze(mean(data_allBase_db(mu_f_bound(1):mu_f_bound(2), :)));
            mu_allBase_mean_db_stance1        = mu_allBase_mean_db(t_bound(1):t_bound(2)-1);
            mu_allBase_mean_db_leftMidSwing     = mu_allBase_mean_db(floor((t_bound(2) + t_bound(3))/2):t_bound(3)); %starts at mid swing to HS
            mu_allBase_mean_db_stance2        = mu_allBase_mean_db(t_bound(3):t_bound(4)-1);
            mu_allBase_mean_db_rightMidSwing    = mu_allBase_mean_db(floor((t_bound(4) + t_bound(5))/2):t_bound(5)); %starts at mid swing to HS

            % Compute allBase means for percent data
            beta_allBase_mean_pChange               = squeeze(mean(data_allBase_pChange(beta_f_bound(1):beta_f_bound(2), :)));
            beta_allBase_mean_pChange_left_swing    = beta_allBase_mean_pChange(t_bound(2):t_bound(3));
            beta_allBase_mean_pChange_leftMid2Stance  = beta_allBase_mean_pChange(floor((t_bound(2) + t_bound(3))/2):t_bound(4)-1);
            beta_allBase_mean_pChange_right_swing   = beta_allBase_mean_pChange(t_bound(4):t_bound(5));
            beta_allBase_mean_pChange_rightMid2Stance  = [ beta_allBase_mean_pChange(floor((t_bound(4) + t_bound(5))/2):t_bound(5)) ...
                beta_allBase_mean_pChange(t_bound(1):t_bound(2)-1)];

            mu_allBase_mean_pChange                = squeeze(mean(data_allBase_pChange(mu_f_bound(1):mu_f_bound(2), :)));
            mu_allBase_mean_pChange_stance1        = mu_allBase_mean_pChange(t_bound(1):t_bound(2)-1);
            mu_allBase_mean_pChange_leftMidSwing     = mu_allBase_mean_pChange(floor((t_bound(2) + t_bound(3))/2):t_bound(3)); %starts at mid swing to HS
            mu_allBase_mean_pChange_stance2        = mu_allBase_mean_pChange(t_bound(3):t_bound(4)-1);
            mu_allBase_mean_pChange_rightMidSwing    = mu_allBase_mean_pChange(floor((t_bound(4) + t_bound(5))/2):t_bound(5)); %starts at mid swing to HS

            % Compute selfBase means for dB data
            beta_selfBase_mean_db                = squeeze(mean(data_selfBase_db(beta_f_bound(1):beta_f_bound(2), :)));
            beta_selfBase_mean_db_left_swing     = beta_selfBase_mean_db(t_bound(2):t_bound(3));
            beta_selfBase_mean_db_leftMid2Stance  = beta_selfBase_mean_db(floor((t_bound(2) + t_bound(3))/2):t_bound(4)-1);
            beta_selfBase_mean_db_right_swing    = beta_selfBase_mean_db(t_bound(4):t_bound(5));
            beta_selfBase_mean_db_rightMid2Stance  = [ beta_selfBase_mean_db(floor((t_bound(4) + t_bound(5))/2):t_bound(5)) ...
                beta_selfBase_mean_db(t_bound(1):t_bound(2)-1)];

            mu_selfBase_mean_db                = squeeze(mean(data_selfBase_db(mu_f_bound(1):mu_f_bound(2), :)));
            mu_selfBase_mean_db_stance1        = mu_selfBase_mean_db(t_bound(1):t_bound(2)-1);
            mu_selfBase_mean_db_leftMidSwing     = mu_selfBase_mean_db(floor((t_bound(2) + t_bound(3))/2):t_bound(3));
            mu_selfBase_mean_db_stance2        = mu_selfBase_mean_db(t_bound(3):t_bound(4)-1);
            mu_selfBase_mean_db_rightMidSwing    = mu_selfBase_mean_db(floor((t_bound(4) + t_bound(5))/2):t_bound(5));

            % Compute selfBase means for percent data
            beta_selfBase_mean_pChange               = squeeze(mean(data_selfBase_pChange(beta_f_bound(1):beta_f_bound(2), :)));
            beta_selfBase_mean_pChange_left_swing    = beta_selfBase_mean_pChange(t_bound(2):t_bound(3));
            beta_selfBase_mean_pChange_leftMid2Stance  = beta_selfBase_mean_pChange(floor((t_bound(2) + t_bound(3))/2):t_bound(4)-1);
            beta_selfBase_mean_pChange_right_swing   = beta_selfBase_mean_pChange(t_bound(4):t_bound(5));
            beta_selfBase_mean_pChange_rightMid2Stance  = [beta_selfBase_mean_pChange(floor((t_bound(4) + t_bound(5))/2):t_bound(5)) ...
                beta_selfBase_mean_pChange(t_bound(1):t_bound(2)-1)];

            mu_selfBase_mean_pChange                = squeeze(mean(data_selfBase_pChange(mu_f_bound(1):mu_f_bound(2), :)));
            mu_selfBase_mean_pChange_stance1        = mu_selfBase_mean_pChange(t_bound(1):t_bound(2)-1);
            mu_selfBase_mean_pChange_leftMidSwing   = mu_selfBase_mean_pChange(floor((t_bound(2) + t_bound(3))/2):t_bound(3));
            mu_selfBase_mean_pChange_stance2        = mu_selfBase_mean_pChange(t_bound(3):t_bound(4)-1);
            mu_selfBase_mean_pChange_rightMidSwing  = mu_selfBase_mean_pChange(floor((t_bound(4) + t_bound(5))/2):t_bound(5));

            %%% set up table
            counter = counter + 1;
            T_tmp = table();
            T_tmp.subject       = {curr_subj};
            T_tmp.group         = {tfdata.group};
            T_tmp.subgroup      = tfdata.subgroup;
            T_tmp.cue           = {tfdata.labels{c}};
            T_tmp.age           = tfdata.age;
            T_tmp.weight        = tfdata.weight;
            T_tmp.height        = tfdata.height;
            T_tmp.moca          = tfdata.MOCA;
            T_tmp.updrs3        = tfdata.UPDRS3;
            T_tmp.pas           = tfdata.PAS;

            switch tfdata.roi

                case 'LSM'
                    %%% allBase
                    T_tmp.betaSupp_allBase_pChange      = mean(beta_allBase_mean_pChange_right_swing);
                    T_tmp.betaSupp_allBase_db           = mean(beta_allBase_mean_db_right_swing);
                    T_tmp.betaMod_allBase_pChange       = mean(beta_allBase_mean_pChange_leftMid2Stance) - mean(beta_allBase_mean_pChange_right_swing);
                    T_tmp.muMod_allBase_pChange         = mean(mu_allBase_mean_pChange_stance2) - mean(mu_allBase_mean_pChange_rightMidSwing);

                    %%% selfBase
                    T_tmp.betaSupp_selfBase_pChange      = mean(beta_selfBase_mean_pChange_right_swing);
                    T_tmp.betaSupp_selfBase_db           = mean(beta_selfBase_mean_db_right_swing);
                    T_tmp.betaMod_selfBase_pChange       = mean(beta_selfBase_mean_pChange_leftMid2Stance) - mean(beta_selfBase_mean_pChange_right_swing);
                    T_tmp.muMod_selfBase_pChange         = mean(mu_selfBase_mean_pChange_stance2) - mean(mu_selfBase_mean_pChange_rightMidSwing);

                    if ~isempty(p_mask)
                        T_tmp.betaSupp_allBase_pChange_masked        = mean(beta_allBase_mean_pChange_masked_right_swing,'omitmissing');
                        T_tmp.betaSupp_selfBase_pChange_masked        = mean(beta_selfBase_mean_pChange_masked_right_swing,'omitmissing');
                    end

                case 'RSM'
                    %%% allBase
                    T_tmp.betaSupp_allBase_pChange      = mean(beta_allBase_mean_pChange_left_swing);
                    T_tmp.betaSupp_allBase_db           = mean(beta_allBase_mean_db_left_swing);
                    T_tmp.betaMod_allBase_pChange       = mean(beta_allBase_mean_pChange_rightMid2Stance) - mean(beta_allBase_mean_pChange_left_swing);
                    T_tmp.muMod_allBase_pChange         = mean(mu_allBase_mean_pChange_stance1) - mean(mu_allBase_mean_pChange_leftMidSwing);
                    
                    %%% selfBase
                    T_tmp.betaSupp_selfBase_pChange     = mean(beta_selfBase_mean_pChange_left_swing);
                    T_tmp.betaSupp_selfBase_db          = mean(beta_selfBase_mean_db_left_swing);
                    T_tmp.betaMod_selfBase_pChange      = mean(beta_selfBase_mean_pChange_rightMid2Stance) - mean(beta_selfBase_mean_pChange_left_swing);
                    T_tmp.muMod_selfBase_pChange         = mean(mu_selfBase_mean_pChange_stance1) - mean(mu_selfBase_mean_pChange_leftMidSwing);

                    if ~isempty(p_mask)
                        T_tmp.betaSupp_allBase_pChange_masked = mean(beta_allBase_mean_pChange_masked_left_swing,'omitmissing');
                        T_tmp.betaSupp_selfBase_pChange_masked = mean(beta_selfBase_mean_pChange_masked_left_swing,'omitmissing');
                    end
            end

            T_tmp.muRMS_allBase_pChange        = rms(mu_allBase_mean_pChange);
            T_tmp.betaRMS_allBase_pChange        = rms(beta_allBase_mean_pChange);
            T_tmp.muRMS_selfBase_pChange        = rms(mu_selfBase_mean_pChange);
            T_tmp.betaRMS_selfBase_pChange        = rms(beta_selfBase_mean_pChange);

            mean_data     = mean(tfdata.gait_table{c}(:,3:end),'omitnan');
            mean_data.Properties.VariableNames = strcat("mean_", mean_data.Properties.VariableNames);

            cv_data     = std(tfdata.gait_table{c}(:,3:end),0,'omitnan')  ./ mean(tfdata.gait_table{c}(:,3:end),'omitnan');
            cv_data.Properties.VariableNames = strcat("cv_", cv_data.Properties.VariableNames);

            T(counter,:) = [T_tmp mean_data cv_data];
        end
    end

    tableName = fullfile(saveDir,[clusterLabel{clu} '_erspData.csv']);
    if exist(tableName, 'file') == 2
        delete(tableName); % Delete the file
        disp('Overwrite previous csv');
    end
    writetable(T,tableName);
    disp('done.')
end