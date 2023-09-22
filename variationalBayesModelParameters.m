function [popmean, pcs_psi, stddev_psi, pcs_lambda, stddev_lambda] = variationalBayesModelParameters(patientData, nu, alpha, npcs_intrapatient, npcs_interpatient)
%VARIATIONALBAYESMODELPARAMETERS Calculate parameters/perform training for
%the variational Bayes organ deformation model

% function [popmean, pcs_intra, stddev_intra, pcs_inter, stddev_inter] = variationalBayesModelParameters(patientData, npcs_intrapatient, npcs_interpatient)
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
% nu:           "Degrees of freedom" parameter for the inverse Wishart
%               distribution
% alpha (optional): Multiplier for the inter-patient uncertaintly Lambda,
%                   see reference for details.
%
%  npcs_intrapatient (optional): Number of intra-patient principal components to
%               output. If not given, all PCs with nonzero variance are
%               given. This is equal to the total number of shapes minus
%               the number of patients.
%  npcs_interpatient (optional): Number of inter-patient principal components to
%               output. If not given, all PCs with nonzero variance are
%               given. 
%  
%
% Output arguments:
%  popmean:       Population mean shape vector
%  pcs_psi:       Principal components of Psi, related to the intra-patient
%                 covariance matrix. Each column represents one principal component
%  stddev_psi:    Vector containing the standard deviation of each
%                 principal component of Psi. Ordered in descending 
%                 fashion.
%  pcs_lambda:     Matrix of inter-patient principal components. Each column 
%                 represents one principal component.
%  stddev_lambda:  Vector containing the standard deviation of each
%                 inter-patient principal component. Ordered in descending 
%                 fashion.

if nargin < 3 || isempty(alpha)
    alpha = 1;
end
if nargin < 4
    npcs_intrapatient = [];
end
    
%Intra patient modes
[popmean, pcs_psi, stddev_psi, individualMeans] = populationModelParameters(patientData, npcs_intrapatient);

%The rest is the calculation of the bias corrected Inter-patient PCA based
%on appendix E in Rørtveit et al., 2023
M = size(individualMeans, 2);
%Inter-patient PCA, not bias corrected
[U, S] = svd(1/sqrt(M-1)*(individualMeans-popmean), 0);

%Intra/Inter data matrices
DL = U*S;
DP = pcs_psi*diag(stddev_psi);

%Normalization factor C
C = 0;
for i = 1:M
    C = C + 1/length(patientData(i).contourPoints);
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

pcs_lambda = W;
if nargin < 4
    npcs_interpatient = nnz(L > 0);
end

pcs_lambda = pcs_lambda(:, 1:npcs_interpatient);
stddev_lambda = sqrt(L(1:npcs_interpatient));
   
    

if any(stddev_lambda < 0)
    error('Negative eigenvalues found. Data matrix is complex');
end

stddev_psi = stddev_psi*sqrt(nu);
stddev_lambda = stddev_lambda*sqrt(alpha);
end
