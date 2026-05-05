function EEG  = extract_straightSegs(EEG)
% remove segments between walking segments that are not in standing

% rhs1 = round([EEG.event(contains({EEG.event.type},'HS1')).latency]);
% rhs2 = round([EEG.event(contains({EEG.event.type},'HS2')).latency]);
% rhs = sort(union(rhs1, rhs2));

rhs = round([EEG.event(contains({EEG.event.type},'HS')).latency]);
keep = zeros(size(EEG.times));
turns = round([EEG.event(contains({EEG.event.type},'turnAPDM')).latency]);
for iInt=1:length(rhs)-1
    current_array = rhs(iInt):rhs(iInt+1);
    if ~any(ismember(turns,current_array))
        keep(current_array) = 1;
    else
        % keep 1s at the beginning and end of the turn. 
        matchIdx = find(ismember(turns, current_array));
        cutoff = turns(matchIdx(end));
        % beginning 1st
        keep(current_array(1:EEG.srate)) = 1;
        % end 1st
        valuesAfterCutoff = current_array(current_array > cutoff);
        if length(valuesAfterCutoff) >= EEG.srate
            valuesAfterCutoff = valuesAfterCutoff(end-EEG.srate+1:end);
            keep(valuesAfterCutoff) = 1;
        end
    end
end

% keep standing
stand_start = [EEG.event(find(contains({EEG.event.type},'standing_start'))).latency];
stand_end =  [EEG.event(find(contains({EEG.event.type},'standing_end'))).latency];
for iInt=1:length(stand_start)
    keep(stand_start(iInt):stand_end(iInt)) = 1;
end

% get intervals
keep_intervals = logical2interval(logical(keep));

% add spacing after and before the turn segments cuts so that the
% boundaries don't fall on the HS event
% cushion = round(EEG.srate*.1);
% 
% if keep_intervals(1,1) <= cushion
%     keep_intervals(2:end,1) = keep_intervals(2:end,1) - cushion-1;
% else
%     keep_intervals(:,1) = keep_intervals(:,1) - cushion-1;
% end
% 
% if keep_intervals(end,2) > EEG.pnts - cushion+1
%     keep_intervals(1:end-1,2) = keep_intervals(1:end-1,2) + cushion-1;
% else
%     keep_intervals(1:end,2) = keep_intervals(1:end,2) + cushion-1;
% end

EEG = pop_select(EEG, 'point', keep_intervals);
EEG.etc.straightWalk_sample_mask = keep;
EEG.setname = [EEG.setname '_straightSegs'];

end