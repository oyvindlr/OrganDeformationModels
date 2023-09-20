
function [totmean, patMeans] = computePatientMeans(patients, maxfracs)
%COMPUTEPATIENTMEANS Compute the mean shape vector for each patient, and
%the population mean (i.e. the mean of the individual means)

if nargin < 2
    maxfracs = Inf;
end
p = size(patients(1).contourPoints{1},1) * 3;
patMeans = zeros(p, length(patients));
for i = 1:length(patients)
    patmean = zeros(p,1);
    nFracs = min(length(patients(i).contourPoints), maxfracs);
    for j = 1:nFracs
        patmean = patmean + vec(patients(i).contourPoints{j});
    end
    patmean = patmean/nFracs;
    patMeans(:, i) = patmean;
end
totmean = mean(patMeans, 2);