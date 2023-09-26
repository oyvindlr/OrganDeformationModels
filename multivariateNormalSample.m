function s = multivariateNormalSample(mu, pcs, stddevs, n)
%MULTIVARIATENORMALSAMPLES Draw a sample from a multivariate normal
% distribution.
%
% The covariance matrix is given by a set of principal components and their
% standard deviations.
%
% Input arguments:
%   mu: Mean vector
%   pcs: Matrix of principal components. Each column represents one
%        principal component.
%   stddevs: Vector of standard deviations corresponding to each principal
%            component
%   n (optional): Sample size, i.e. number of realizations to produce. One 
%                 realization is produced if n is not given.
%
% Output arguments:
%   s: Vector/matrix of realizations. Each column represents one
%   realization

if nargin < 4 || isempty(n)
    n = 1;
end

s = mu + pcs*diag(stddevs)*randn(length(stddevs), n);

