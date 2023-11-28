function [patients, shifts] = rigidshift(patients, percent, first, planScan, trainmean)
[totmean, patMeans] = computePatientMeans(patients);
if nargin < 5
    trainmean = totmean;
end
trainmean = vec(trainmean);
[~, ind] = sort(trainmean(:,3));
if nargin < 2
    percent = 100;
end
if nargin < 3
    first = false;
end
if nargin < 4
    planScan = 1;
end

ind = ind(1:round(percent*size(trainmean, 1)/100));

shifts = zeros(length(patients), 3);
for i = 1:length(patients)
    
    if first
        iMean = patients(i).contourPoints{planScan};
    else
        iMean = vec(patMeans(:,i));
    end
    iMean = mean(iMean(ind,:), 1);
    for j = 1:length(patients(i).contourPoints)
        patients(i).contourPoints{j} = patients(i).contourPoints{j} - iMean;
    end
    shifts(i, :) = iMean;
end



