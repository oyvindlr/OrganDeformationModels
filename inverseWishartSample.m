function [pcs_sample, stddev_sample] = inverseWishartSample(psi_pcs, psi_stddev, nu)
%INVERSEWISHARTSAMPLE Draw a single realization of an covariance matrix
%from an inverse Wishart distribution. Note that nu needs to be a natural
%number greater than p+1, where p is the number of principal components of
%psi
%
% Uses the pseudo-inverse for both inversion and "de-inversion".
%
% Input arguments:
%   psi_pcs: Principal components of the scale matrix Psi
%   psi_stddev: Standard deviations of the given principal components
%   nu: Number of degrees of freedom. Must be greater than
%   length(psi_stddev)+1
%   cutoff (optional): 
% Output parameters:
%   pcs_sample: Principal components of the matrix realization
%   stddev_sample: Standard deviation of the principal components.


psi_stddev = psi_stddev(psi_stddev > 1e-9);
p = length(psi_stddev);

if nu <= p+1
    error('inverseWishartSample: Parameter nu must be greater than the number of principal components plus one');
end

psi_stddev = sqrt(nu-p-1)*psi_stddev(1:p);

invD = psi_pcs(:, 1:p)*diag(1./psi_stddev);

%Produce a Wishart-sample using the pseudo-inverse of the covariance matrix
sample = invD*randn(length(psi_stddev), nu);

%PCA
[pcs_sample, stddev_inv] = svd(sample, 0);

%Remove zeroes
stddev_inv = diag(stddev_inv);
stddev_inv = stddev_inv(1:p);
pcs_sample = pcs_sample(:, 1:p);

%de-invert
stddev_sample = 1./stddev_inv;

%Sort PCs by descending order
[stddev_sample, ind] = sort(stddev_sample, 'desc');
pcs_sample = pcs_sample(:, ind);

