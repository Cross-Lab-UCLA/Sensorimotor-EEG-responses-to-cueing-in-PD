%% Function to get time warp values from EEG

function [warpVals, gait_table, gait_events_idx, no_gc_idx ] = getTimeWarp(tmpEEG)
% set up for gait events extraction
L_TO        = [];
L_HS        = [];
R_TO        = [];
R_HS2       = [];
no_gc_idx   = [];
warpVals    = [];

gait_events_idx = nan(length(tmpEEG.epoch),5);
gait_table = table;

for m = 1:length(tmpEEG.epoch)
    L_TO_1 = []; L_HS_1 = []; R_TO_1 = []; R_HS2_1 = [];

    % HS1_idx = find(contains(tmpEEG.epoch(m).eventtype, 'R_HS1'));
    % cond = strsplit(tmpEEG.epoch(m).eventtype{HS1_idx(1)},'_');
    % cond = cond(end);

    L_TO    = cell2mat(tmpEEG.epoch(m).eventlatency(contains(tmpEEG.epoch(m).eventtype, 'L_TO')));
    L_HS    = cell2mat(tmpEEG.epoch(m).eventlatency(contains(tmpEEG.epoch(m).eventtype, 'L_HS')));
    R_TO    = cell2mat(tmpEEG.epoch(m).eventlatency(contains(tmpEEG.epoch(m).eventtype, 'R_TO')));
    R_HS2   = cell2mat(tmpEEG.epoch(m).eventlatency(contains(tmpEEG.epoch(m).eventtype, 'R_HS2')));
    R_HS2_label = tmpEEG.epoch(m).eventtype(contains(tmpEEG.epoch(m).eventtype, 'R_HS2'));
    g_table = tmpEEG.epoch(m).eventgait(contains(tmpEEG.epoch(m).eventtype, 'R_HS2'));

    R_HS2_1   = R_HS2(find(R_HS2 > 0, 1, 'first'));
    R_HS2_1_label   = R_HS2_label(find(R_HS2 > 0, 1, 'first'));

    if ~isempty(R_HS2_1)
        L_TO_1    = L_TO(find(L_TO > 0 & L_TO < R_HS2_1, 1, 'first'));
        L_HS_1    = L_HS(find(L_HS > 0 & L_HS < R_HS2_1, 1, 'first'));
        R_TO_1    = R_TO(find(R_TO > 0 & R_TO < R_HS2_1, 1, 'first'));
    end

    if  isempty(L_TO_1) || isempty(L_HS_1) || isempty(R_TO_1) || isempty(R_HS2_1) || ...
            sum(R_TO(R_TO>0)<R_HS2_1) > 1 || R_HS2_1 < R_TO_1 || L_HS_1 < L_TO_1
        no_gc_idx = [no_gc_idx; m];  % track gait cycles that are invalid.
    else
        gait_events_idx(m, :) = [0 L_TO_1 L_HS_1 R_TO_1 R_HS2_1]; % R_HS1 is ALWAYS at zero latency
        T = g_table{R_HS2 == R_HS2_1}(1,1:end);
        T.Condition = extractAfter(R_HS2_1_label,'HS2_');
        gait_table         = [gait_table;T];
    end

end

gait_events_idx(no_gc_idx, :)   = [];

% if height(gait_table) ~= length(gait_events_idx(:,1))
%     % check if table and gait indexes match up
%     keyboard
% end

warpVals = median(cat(1,gait_events_idx));
end