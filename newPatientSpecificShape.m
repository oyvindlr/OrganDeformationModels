function s = newPatientSpecificShape(mu, pcs, stddevs)
%NEWPATIENTSPECIFICSHAPE Generates a Gaussian random vector with mean mu
% and PCA of the covariance matrix given by pcs and stddevs

s = mu + pcs*(randn(size(stddevs, 2), 1).*stddevs(:));