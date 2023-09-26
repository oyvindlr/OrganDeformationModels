function [mu, pcs, stddevs] = variationalBayesSample(mu0, pcs_psi, stddev_psi, pcs_lambda, stddev_lambda, nu)
%VARIATIONALBAYESSAMPLE Draw a new realization of a mean vector and
%covariance matrix from the distribution used for the variational Bayes
%model 
%
% The mean is drawn from a multivariate Gaussian distribution, the
% covariance matrix is drawn from an inverse Wishart distribution, and they
% are independent.
%
% Input arguments: 
%  mu_0: Mean vector
%  pcs_psi: Principal component matrix of the scale matrix Psi of the
%  inverse Wishart distribution
%  stddev_ps: Standard deviations of the PCs of Psi
%  pcs_lambda: Principal components of the normal distribution
%  stddev_lambda: Standard deviations of the PCs of the normal distribution
%  nu: Number of degrees of freedom for the inverse Wishart distribution
%
% Output arguments:
%  mu: Drawn mean vector
%  pcs: Principal components of the drawn covariance matrix
%  stddevs: Standard deviations of the drawn PCs

mu = multivariateNormalSample(mu0, pcs_lambda, stddev_lambda);

[pcs, stddevs] = inverseWishartSample(pcs_psi, stddev_psi, nu);

