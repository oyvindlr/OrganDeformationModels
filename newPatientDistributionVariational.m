function [mu, pcs, stddevs] = newPatientDistributionVariational(mu0, pcs_psi, stddev_psi, pcs_lambda, stddev_lambda, nu)

D_lambda = pcs_lambda*diag(stddev_lambda);

mu = mu0 + D_lambda*randn(size(D_lambda, 2), 1);

[pcs, stddevs] = inverseWishartSample(pcs_psi, stddev_psi, nu);

