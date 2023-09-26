function [pcs_sample, stddev_sample] = inverseWishartSample(psi_pcs, psi_stddev, nu, cutoff)
%INVERSEWISHARTSAMPLE Draw a single realization of an covariance matrix
%from an inverse Wishart distribution
%
% Uses the pseudo-inverse for both inversion and "de-inversion".
% Note: The degenerate forward Wishart distribution, i.e. with nu < P
%   where PxP is the matrix dimension, will sometimes produce very small
%   eigenvalues. These in turn become very large eigenvalues of the
%   pseudo-inverse Wishart distribution. For this reason we introduce an
%   optional "cutoff" whereby all realizations with eigenvalues larger than
%   cutoff*max(psi_stddev)/sqrt(nu) will be discarded. Typical values for
%   the cutoff would be around 5, or they can be measured from training 
%   data. 
%
% Input arguments:
%   psi_pcs: Principal components of the scale matrix Psi
%   psi_stddev: Standard deviations of the given principal components
%   nu: Number of degrees of freedom
%   cutoff (optional): 
% Output parameters:
%   pcs_sample: Principal components of the matrix realization
%   stddev_sample: Standard deviation of the principal components.
%   cutoff (optional): See description above



invD = psi_pcs*diag(1./psi_stddev);

%Produce a Wishart-sample using the pseudo-inverse of the covariance matrix
sample = invD*randn(length(psi_stddev), nu);

%PCA
[pcs_sample, stddev_inv] = svd(sample, 0);

%de-invert
stddev_sample = 1./diag(stddev_inv);

%remove infinite-valued variances
stddev_sample(size(psi_pcs, 2)+1:end) = 0;

%Sort PCs by descending order
[stddev_sample, ind] = sort(stddev_sample, 'desc');
pcs_sample = pcs_sample(:, ind);

if nargin > 3 && ~isempty(cutoff)
    maxx = cutoff*max(psi_stddev)/sqrt(nu);
    if stddev_sample(1) > maxx
        [pcs_sample, stddev_sample] = inverseWishartSample(psi_pcs, psi_stddev, nu, cutoff);
    end
end
