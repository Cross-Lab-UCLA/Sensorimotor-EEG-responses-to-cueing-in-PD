function [EEG, T] = removedCyclesMoreAffected(EEG)
% Find epoch data that deviates from the median and remove them
%
% T = tracking table
%%

HS1_idx = find(~cellfun('isempty', regexp({EEG.event.type}, '^MoreAffected_.*_HS1.*$')));
rhs1 = round([EEG.event(HS1_idx).latency]);
rhs = unique(rhs1);

% epoch gait cycle data
epoch_points = [];
gfp = [];

for iInt=1:length(rhs)-1
    current_array = rhs(iInt):rhs(iInt+1);
    epoch_points(iInt,:) = [rhs(iInt) rhs(iInt+1)];
    gfp{iInt} = rms(EEG.data(:,current_array)); % global field power
end

% find any segments that contain data point that exceed 3 std from median
data = horzcat(gfp{:});
global_threshold1 = mean(data) + std(data)*3;
global_threshold2 = mean(data) + std(data)*5;
threshold_percent = 0.05;
std_mask1 = []; std_mask2 = [];

for hs = 1:length(gfp)
    num_timepoints  = size(gfp{hs}, 2);
    exceed_mask1     = gfp{hs} > global_threshold1;
    exceed_count1    = sum(exceed_mask1); 
    std_mask1(hs, :) = exceed_count1 > (threshold_percent * num_timepoints); % Flag if >0.05% exceed

    exceed_mask2     = gfp{hs} > global_threshold2;
    exceed_count2    = sum(exceed_mask2); 
    std_mask2(hs, :) = exceed_count2 > 0; % Flag if any point is over 5 std of median 
end
std_mask_combined = std_mask1|std_mask2;
rm_idx = [];
rm_idx = find(std_mask_combined ==1)';

total_cycles = length(gfp); 
if isempty(rm_idx)
    T.subject = EEG.subject;
    T.total_cycles = total_cycles;
    T.outlier_cycles_flagged = 0;
    T.cycles_flag_percent = 0;
    T.total_pnts = length(EEG.times);
    T.outlier_pnts_removed = 0;
    T.pnt_removed_percent = 0;
else
    disp(['Flagging ' num2str(length(rm_idx)) ' gait cycles for removal']);
    T.subject = EEG.subject;
    T.total_cycles = total_cycles; % this number will change after epoching, gait cycles with boundary between will be removed
    T.outlier_cycles_flagged = length(rm_idx);
    T.cycles_flag_percent = T.outlier_cycles_flagged/total_cycles * 100;

    % get removal array
    rm_array = epoch_points(rm_idx,:);
    % add cushion so that boundary doesn't land on the HS index
    rm_array(:,1) = rm_array(:,1) + 1; 
    rm_array(:,2) = rm_array(:,2) - 1; 

    % track points removed
    T.total_pnts = length(EEG.times);
    T.outlier_pnts_removed = sum(diff(rm_array')) + length(diff(rm_array'));
    T.pnt_removed_percent = T.outlier_pnts_removed/length(EEG.times) * 100;

    % remove segments for EEG
    EEG = pop_select(EEG, 'rmpoint', rm_array);
    EEG = eeg_checkset(EEG);
end
end