function intervals = logical2interval(vector, minlength)
% turns logical timeseries vector of length n to a 2xn matrix of onset (row 1) 
% and offset (row 2) samples of true segments that are longer than minlength
% minlength and intervals are in samples

intervals = reshape(find(diff([false vector false])),2,[])';
intervals(:,2) = intervals(:,2)-1;

% remove small intervals
if exist('minlength','var') && isscalar(minlength) && ~isempty(intervals)
    smallIntervals = diff(intervals')' < minlength;
    for iInterval = find(smallIntervals)'
        sample_mask(intervals(iInterval,1):intervals(iInterval,2)) = 0;
    end
    intervals(smallIntervals,:) = [];
end
