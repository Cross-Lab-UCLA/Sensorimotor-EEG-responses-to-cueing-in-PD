%% create csv from Targeted Clusters
%
% Run after A9. Extract more- and less- affect side
% variables for statistical analysis.
%
% LM 072525
%%

clear all; clc; close all;
if ispc
    mainDir     = 'C:\Git\DoD-Gait';
else
    mainDir = '/Users/Leo/Git/DoD-Gait';
end
dataDir = fullfile(mainDir,'data');
saveDir = fullfile(mainDir,'results','processed_data');

% get tables
path = fullfile(mainDir,'results','processed_data');

RSM = readtable(fullfile(path,['RSM_erspData.csv']));
LSM = readtable(fullfile(path,['LSM_erspData.csv']));
affectedside = readtable(fullfile(dataDir,'affected_side.csv'));

T_more = table();
T_less = table();

%Note: the more affected side is relative to the body and thus
%contralateral to the brain

for sub = 1:height(affectedside)

    subject = affectedside.subject{sub};
    moreAffectedSide = affectedside.affected_side_label{sub};

	disp(['Running subject: ' subject]);
	if strcmpi(moreAffectedSide,'R')
        disp('...marking LSM as the MORE affected side');
		% get LSM and assign R gait paramters as more affected side
        rows = LSM(strcmp(LSM.subject, subject),:);
        if isempty(rows)
            disp('...LSM dipole not found...skipping');
        else
            varNames = rows.Properties.VariableNames;
            idx = endsWith(varNames, '_R');
            varNames(idx) = strrep(varNames(idx), '_R', '_moreAffected');
            rows.Properties.VariableNames = varNames;
            idx_Remove = endsWith(varNames, '_L');
            rows(:, idx_Remove) = [];
            rows.cluster = repmat('LSM', height(rows), 1);
            rows.body_side = repmat('R', height(rows), 1);
            T_more = [T_more; rows];
        end

        % get RSM and assign L gait paramters as less affected side
        disp('...marking RSM as the less affected side');
		rows = RSM(strcmp(RSM.subject, subject),:);
        if isempty(rows)
            disp('...RSM dipole not found...skipping');
        else
            varNames = rows.Properties.VariableNames;
            idx = endsWith(varNames, '_L');
            varNames(idx) = strrep(varNames(idx), '_L', '_lessAffected');
            rows.Properties.VariableNames = varNames;
            idx_Remove = endsWith(varNames, '_R');
            rows(:, idx_Remove) = [];
            rows.cluster = repmat('RSM', height(rows), 1);
            rows.body_side = repmat('L', height(rows), 1);
            T_less = [T_less; rows];
        end

    elseif strcmpi(moreAffectedSide,'L')
		disp('...marking RSM as the MORE affected side');

        % get RSM and assign L gait paramters as more affected side
        rows = RSM(strcmp(RSM.subject, subject),:);
        if isempty(rows)
            disp('...RSM dipole not found...skipping');
        else
            varNames = rows.Properties.VariableNames;
            idx = endsWith(varNames, '_L');
            varNames(idx) = strrep(varNames(idx), '_L', '_moreAffected');
            rows.Properties.VariableNames = varNames;
            idx_Remove = endsWith(varNames, '_R');
            rows(:, idx_Remove) = [];
            rows.cluster = repmat('RSM', height(rows), 1);
            rows.body_side = repmat('L', height(rows), 1);
            T_more = [T_more; rows];
        end

        % get LSM and assign R gait paramters as less affected side
        disp('...marking LSM as the less affected side');
        rows = LSM(strcmp(LSM.subject, subject),:);
        if isempty(rows)
            disp('...LSM dipole not found...skipping');
        else
            varNames = rows.Properties.VariableNames;
            idx = endsWith(varNames, '_R');
            varNames(idx) = strrep(varNames(idx), '_R', '_lessAffected');
            rows.Properties.VariableNames = varNames;
            idx_Remove = endsWith(varNames, '_L');
            rows(:, idx_Remove) = [];
            rows.cluster = repmat('LSM', height(rows), 1);
            rows.body_side = repmat('R', height(rows), 1);
            T_less = [T_less; rows];
        end
    
    else
        keyboard
    end
end

% save more affected
tableName = fullfile(saveDir,['moreAffected_erspData.csv']);
if exist(tableName, 'file') == 2
    delete(tableName); % Delete the file
    disp('Overwrite previous csv');
end
writetable(T_more,tableName);

% save less affected
tableName = fullfile(saveDir,['lessAffected_erspData.csv']);
if exist(tableName, 'file') == 2
    delete(tableName); % Delete the file
    disp('Overwrite previous csv');
end
writetable(T_less,tableName);

disp('done');