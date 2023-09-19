% function [pcs, stddev, avg] = pca(dataMatrix, npcs)
%
% Principal component analysis of the sample represented by the columns of
% the data matrix.
% Input:
%  dataMatrix: Matrix containing the data, each column is data
%              point/realization
%  npcs:       Number of principal components to otput (optional). If not
%              given, the number of pcs is sample size minus one.
% Output:
%  pcs:        Prinicipal components matrix. Each column represents one principal
%              component.
%  stddev:     array of standard deviations of each principal component. In
%              descending order.
%  avg:        PCA mean vector, i.e. mean of the columns of dataMatrix


function [pcs, stddev, avg] = pca(dataMatrix, npcs)

m = size(dataMatrix, 2);
if nargin < 2
    npcs = m - 1;
end

avg = mean(dataMatrix, 2);
[pcs, S] = svd(dataMatrix - avg, 0);
stddev = sqrt(diag(S)./(m-1));
pcs = pcs(:, 1:npcs);
stddev = stddev(1:npcs);

