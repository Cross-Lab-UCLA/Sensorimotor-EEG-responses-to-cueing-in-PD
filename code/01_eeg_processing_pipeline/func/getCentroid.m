function [centroid, mean_dist,std_dist,outlier_idx] = getCentroid(XYZ)
% get centriod info from a x by 3 xyz matrix
    centroid = mean(XYZ, 1);
    distances = pdist2(XYZ,centroid);
    mean_dist = mean(distances);
    std_dist = std(distances);
    threshold = mean_dist + 3 * std_dist;
    outlier_idx = find(distances > threshold);
end