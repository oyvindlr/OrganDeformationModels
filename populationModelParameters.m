% function [popmean, pcs, standardDev, individualMeans] = populationModelParameters(patientData, nComponents)
% Calculate the parameters, i.e. the prinicpal components and their
% variance, as well as the population mean, of a population deformation
% model.
% Input arguments:
%  patientData: struct array containing the point data for all patients.
%               each element in the array represents one patient, and each
%               patient must have a set of structures represented by the
%               cell array "contourPoints". In this cell array, each cell
%               represents an organ shape from that patient. The shape is
%               represented by an N by 3 array, where N is the number of
%               points, each of which has three spatial coordinates.
%
%  nComponents: Number of principal components to output. If not given, all
%               PCs with nonzero variance are given. This is equal to 
%               the total number of shapes minus the number of patients.         
%  
%
% Output arguments:
%  popmean: Population mean shape vector
%  pcs: Matrix of principal components. Each column represents one principal
%  component
%  standardDev: vector containing the standard deviation of each principal
%  component. Ordered in descending fashion.
%  individualMeans: Matrix where each row is the mean vector of one patient

function [popmean, pcs, standardDev, individualMeans] = populationModelParameters(patientData, nComponents)



[popmean, individualMeans] = computePatientMeans(patientData);


npatients = length(patientData);
p = size(patientData(1).contourPoints{1},1) * 3;
count = 1;
fullMat = zeros(p,0);

%Gather all contourpoint data into one big matrix
for i = 1:npatients
    pat = patientData(i);
    nFracs = length(pat.contourPoints);
    patMat = zeros(p, nFracs);
    for j = 1:nFracs        
        patMat(:,j) =  pat.contourPoints{j}(:);
    end
    %subtract individual mean from individual data matrix
    patMat = patMat - repmat(individualMeans(:, i), 1, size(patMat, 2));
    fullMat(:, count:count+nFracs-1) = patMat;
    count = count + nFracs;
end

%singular value decomposition to calculate PCA
[pcs, S] = svd(fullMat, 0);

S = diag(S);

if nargin < 2 || isempty(nComponents)
    nComponents = size(fullMat, 1) - npatients;
end
%Remove excess modes
if size(pcs, 2) > nComponents
    pcs = pcs(:, 1:nComponents);
    S = S(1:nComponents);
end
%Calculate standard deviation
standardDev = S./sqrt(size(fullMat,2)-1);