function [popmean, pcs_intra, stddev_intra, pcs_inter, stddev_inter] = variationalBayesModelParameters(patientData, npcs_intrapatient, npcs_interpatient)
%VARIATIONALBAYESMODELPARAMETERS Summary of this function goes here

%  function [popmean, pcs_intra, stddev_intra, pcs_inter, stddev_inter] = variationalBayesModelParameters(patientData, npcs_intrapatient, npcs_interpatient)
% Calculate the parameters, i.e. the intra-patient and inter-patient 
% prinicpal components and theirvariance, as well as the population mean,
% for the variational Bayes deformation model in Rørtveit et al., 2023. 
% This function uses the bias-correction described in that paper to compute
% the inter-patient principal components.
%
% Input arguments:
%  patientData: struct array containing the point data for all patients.
%               each element in the array represents one patient, and each
%               patient must have a set of structures represented by the
%               cell array "contourPoints". In this cell array, each cell
%               represents an organ shape from that patient. The shape is
%               represented by an N by 3 array, where N is the number of
%               points, each of which has three spatial coordinates.
%
%  npcs_intrapatient: Number of intra-patient principal components to
%               output. If not given, all PCs with nonzero variance are
%               given. This is equal to the total number of shapes minus
%               the number of patients.
%  npcs_interpatient: Number of inter-patient principal components to
%               output. If not given, all PCs with nonzero variance are
%               given. 
%  
%
% Output arguments:
%  popmean: Population mean shape vector
%  pcs_intra:     Matrix of intra-patient principal components. Each column 
%                 represents one principal component
%  stddev_intra:  Vector containing the standard deviation of each
%                 intra-patient principal component. Ordered in descending 
%                 fashion.
%  pcs_inter:     Matrix of inter-patient principal components. Each column 
%                 represents one principal component.
%  stddev_inter:  Vector containing the standard deviation of each
%                 inter-patient principal component. Ordered in descending 
%                 fashion.

if nargin < 2
    npcs_intrapatient = [];
end
    
%Intra patient modes
[popmean, pcs_intra, stddev_intra, individualMeans] = populationModelParameters(patientData, npcs_intrapatient);

%The rest is the calculation of the bias corrected Inter-patient PCA based
%on appendix E in Rørtveit et al., 2023
M = size(individualMeans, 2);
%Inter-patient PCA, not bias corrected
[U, S] = svd(1/sqrt(M-1)*(individualMeans-popmean), 0);

%Intra/Inter data matrices
DL = U*S;
DP = pcs_intra*diag(stddev_intra);

%Normalization factor C
C = 0;
for i = 1:M
    C = C + 1/length(patients(i).contourPoints);
end
C = C/M;


A = [DL sqrt(C)*DP];
B = [DL -sqrt(C)*DP];
[W, L] = eig(B'*A);
L = real(diag(L));
[L, ind] = sort(L, 'desc');
W = real(W);
W = W(:, ind);

W = A*W;
for i = 1:size(W, 2)
    W(:, i) = W(:, i)/norm(W(:, i));
end

pcs_inter = W;
stddev_inter = diag(sqrt(L));
if nargin < 2
    npcs_interpatient = nnz(stddev_inter > 0);
end

pcs_inter = pcs_inter(:, 1:npcs_interpatient);
stddev_inter = stddev_inter(1:npcs_interpatient);
   
    

if any(stddev_inter < 0)
    error('Negative eigenvalues found. Data matrix is complex');
end

end
