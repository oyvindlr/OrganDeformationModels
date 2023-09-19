%function [avg, pcs, stddev] = patientSpecificModelParameters(data, npcs)
%
% Calculate the parameters, i.e. the prinicpal components and their
% variance, as well as the patient mean, of a patient-specific deformation
% model.
% Input arguments:
%  data:        cell arraystruct array containing the point data for all
%               shapes for this patient. Each cell represents an organ 
%               shape from the patient. The shape is represented by an
%               N by 3 array, where N is the number of points, each of
%               which has three spatial coordinates.
%
%  nComponents: Number of principal components to output. If not given, all
%               PCs with nonzero variance are given. This is equal to 
%               the total number of shapes minus one.         
%  
%
% Output arguments:
%  avg: Mean shape vector
%  pcs: Matrix of principal components. Each column represents one principal
%  component
%  standardDev: vector containing the standard deviation of each principal
%  component. In descending order.

function [avg, pcs, stddev] = patientSpecificModelParameters(data, npcs)

dataMat = zeros(size(data{1}, 1)*3, size(data));
for i = 1:length(data)
    dataMat(:, i) = data{i}(:);
end
if nargin > 1
    [pcs, stddev, avg] = pca(dataMat, npcs);
else
    [pcs, stddev, avg] = pca(dataMat);
end
    
    