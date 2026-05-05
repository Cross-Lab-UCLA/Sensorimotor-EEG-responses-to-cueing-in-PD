function dataRank = getRank(tmpdata)
% modified from Makoto and pop_runica 
% https://sccn.ucsd.edu/wiki/Makoto%27s_preprocessing_pipeline#What_is_rank.3F_.28Updated_06.2F06.2F2016.29

rank1 = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); 
% this is the input that get passed into the getrank function in pop_runica

covarianceMatrix = cov(tmpdata', 1);
[E, D] = eig(covarianceMatrix);
rankTolerance = 1e-7;
rank2=sum (diag(D) > rankTolerance);

if rank1 ~= rank2
    fprintf('Checking rank computation inconsistency (%d vs %d)\n', rank1, rank2);
    %tmprank2 = max(tmprank, tmprank2);
end

dataRank = min(rank1, rank2);

end