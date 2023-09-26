function [mu, pcs, stddevs] = normalInverseWishartSample(mu_0, kappa, psi_pcs, psi_stddevs, nu)
%NORMALINVERSEWISHARTSAMPLE Get a sample from a normal-inverse-Wishart
%distribution

% Gets a mean vector and covariance matrix drawn from an NIW-distribution.
% The covariance matrix is given by its principal components and their
% standard deviations.

% Input arguments:
%   mu_0: Mean vector 
%   kappa: kappa parameter (sometimes called lambda) of the NIW distribution
%   psi_pcs: Principal components of the scale matrix Psi
%   psi_stddevs: Standard deviations of the principal components
%   nu: Degrees of freedom of the inverse Wishart distribution
%
% Output arguments:
%   mu: Sample of the mean vector
%   pcs: Principal components of the covariance matrix sample
%   stddevs: Standard deviations of the PCA of the covariance matrix

[pcs, stddevs] = inverseWishartSample(psi_pcs, psi_stddevs, nu);

stddevs_mean = sqrt(1/kappa)*stddevs;

mu = multivariateNormalSample(mu_0, pcs, stddevs_mean);