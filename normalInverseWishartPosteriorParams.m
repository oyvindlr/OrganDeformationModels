function [patmean, pcs, stddev, kappa_p, nu_p] = normalInverseWishartPosteriorParams(patientData,popmean, pcs, stddev, kappa, nu)
%NORMALINVERSEWISHARTPOSTERIORPARAMS Get the point estimates for the
%posterior distribution of a normal inverse gamma model
%
%function [patmean, pcs, stddev] = normalInverseWishartPosteriorParams(patientData,popmean, pcs, stddev, kappa, nu)
% Input arguments:
%   patientData: Array of vectorized coordinates of organ shapes for the
%                patient, or cell array of non-vectorized coordinate
%                matrices for the same
%  popmean:      Population mean shape vector
%  pcs:          Population principal components
%  stddev:       Standard deviations of the population principal components
%  kappa:        kappa parameter for the NIW-prior
%  nu:           nu-parameter for the NIW-prior
%
% Output arguments:
%  patmean: Patient specific mean vector (point estimate)
%  pcs:     Patient specific posterior principal components
%  stddev:  Standard deviations of the principal components, in descending
%           order
% kappa_p: posterior value of the kappa parameter
% nu_p: posterior value of the nu parameter

if iscell(patientData)
    patientData = cellToArr(patientData);
end
J = size(patientData, 2);

p = length(stddev);

nu_p = nu + J;

kappa_p = kappa + J;

m = mean(patientData, 2);


patmean = 1/kappa_p*(kappa*popmean + J*m);

D = pcs*diag(stddev);
S = patientData-m;

expansion = [sqrt(kappa*J/kappa_p)*(m-popmean) S];


Dstar = sqrt(1/(nu_p-p-1))*[sqrt(nu-p-1)*D expansion];

[pcs, stddev] = svd(Dstar, 0);

pcs = pcs(:, 1:end-1);

stddev = diag(stddev);
stddev = stddev(1:end-1);

end


function arr = cellToArr(celldata)
arr = zeros(size(celldata{1}, 1)*3, length(celldata));
for i = 1:length(celldata)
    arr(:, i) = celldata{i}(:);
end    
end