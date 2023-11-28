%Create some data to run examples with
[patientData, faces] = createTestPatients;
P = 3; %number of principal components to use

%%
%Alternatively, download actual rectum data from https://doi.org/10.18710/DKVPIJ
% and run the following
load('patientdata.mat');
%Rigid registration between patients using "center of gravity" and only
%using the 50% most caudal points.
percentage_caudal_points = 50;
[patientData, shifts] = rigidshift(patientData, percentage_caudal_points);
P = 12; %number of principal components to use

%%
%Compute the parameters (mean and PCA) for the population model
[popmean, ppcs, pstd, indivmean] = populationModelParameters(patientData);

%Show the mean shape and +/-one standard deviation of the first mode
displayShape(popmean, 'b', 0.4);
hold on;
displayShape(popmean + ppcs(:, 1)*pstd(1), 'g', 0.5);
displayShape(popmean - ppcs(:,1)*pstd(1), 'r', 0.4);
hold off

%%
%Use maximum likelihood to find kappa and nu for the NIW-model using
%posterior predictive sampling

%f = @(x)niw_costfunc(x, patientData, popmean, ppcs(:,1:P), pstd(1:P), 2);
%Alternatively, get the parameters for a point estimate like this:
f = @(x)niw_point_costfunc(x, patientData, popmean, ppcs(:,1:P), pstd(1:P), 1);

%You could also use fminsearch
v = quasinewton(f, [0.2;P+7],'verbose', true, 'linesearchmaxiter', 100);
kappa = v(1);
nu = round(v(2));

%%
%Generate the parameters for a new random patient, then generate some
%samples for that patient using the NIW model
[mu_new, pcs_new, stddevs_new] = normalInverseWishartSample(popmean, kappa, ppcs(:,1:P), pstd(1:P), nu);%Parameters (mean and covariance matrix) for a new random patient

%Generate two samples from the new patient
smp = multivariateNormalSample(mu_new, pcs_new, stddevs_new, 2);

%Display the two samples
hold off
displayShape(smp(:, 1), 'g', 0.7);
hold on
displayShape(smp(:, 2), 'b', 0.7);

%%
%Calculate the posterior parameters for a specific patient
%Generate random data for a given patient, based on the posterior
%distribution, using the NIW-model

% Calculate the posterior parameters given some patient data.
% The posterior parameters can also be used as point estimates for the mean
% and covariance matrix (represented here by principal components) for that patient
patient_no = 3;
[mean_p, pcs_p, stddev_p, kappa_p, nu_p] = normalInverseWishartPosteriorParams(patientData(patient_no).contourPoints(1:2), popmean, ppcs(:, 1:P), pstd(1:P), kappa, nu);
%%
%Generate patient-specific samples using point estimates of the parameters
%(method 1)
%Generate two samples from the new patient
smp = multivariateNormalSample(mean_p, pcs_p, stddev_p, 2);

%Display the two samples
hold off
displayShape(smp(:, 1), 'g', 0.7);
hold on
displayShape(smp(:, 2), 'b', 0.7);
%%
%Generate patient-specific samples using the posterior predictive distribution
%(method 2)
for i = 1:2
    %Get randomly generated parameters from the posterior distribution
    [mu_r, pcs_r, stddevs_r] = normalInverseWishartSample(mean_p, kappa_p, pcs_p, stddev_p, nu_p);
    %Get a multivariate random sample from the random parameters
    smp(:, i) = multivariateNormalSample(mu_r, pcs_r, stddevs_r, 1);
end
%Display the two samples
hold off
displayShape(smp(:, 1), 'g', 0.7);
hold on
displayShape(smp(:, 2), 'b', 0.7);
    

%%
%Calculate the parameters for a variational Bayes model
[popmean, pcs_psi, stddev_psi, pcs_lambda, stddev_lambda] = variationalBayesModelParameters(patientData);

%Show the mean shape and +/-one standard deviation of the first INTER-PATIENT mode
hold off
displayShape(popmean, 'b', 0.4);
hold on;
displayShape(popmean + pcs_lambda(:, 1)*stddev_lambda(1), 'g', 0.5);
displayShape(popmean - pcs_lambda(:,1)*stddev_lambda(1), 'r', 0.4);
hold off
%%
%Create a new patient from the variational model and generate a few random
%samples
[mu_new, pcs_new, stddev_new] = variationalBayesSample(popmean, pcs_psi(:, 1:P), stddev_psi(1:P), pcs_lambda(:, 1:P), stddev_lambda(1:P), nu);

%Generate two samples from the new patient
smp = multivariateNormalSample(mu_new, pcs_new, stddev_new, 2);

%Display the two samples
hold off
displayShape(smp(:, 1), 'g', 0.7);
hold on
displayShape(smp(:, 2), 'b', 0.7);
%%
%Calculate the patient-specific posterior parameters of the variational
%Bayes model. We need to assign a few extra parameter values
delta_psi = 20000;
delta_lambda = 20000;
D_psi = pcs_psi*diag(stddev_psi);
D_lambda = pcs_lambda*diag(stddev_lambda);

patData = [vec(patientData(patient_no).contourPoints{1}) vec(patientData(patient_no).contourPoints{2})];

[imean, psipcs, psistddev, psireg, lambdapcs, lambdastddev, lambdareg] = ...
    fastVariational(popmean, D_psi(:, 1:P),D_lambda(:, 1:P), delta_psi, delta_lambda, patData, nu);

%Display the expected mean shape
displayShape(imean);
%%
%Generate random shapes for the specific patient using point-estimates of
%the parameters

%Generate two samples from the new patient
smp = multivariateNormalSample(imean, psipcs, psistddev, 2);

%Display the two samples
hold off
displayShape(smp(:, 1), 'g', 0.7);
hold on
displayShape(smp(:, 2), 'b', 0.7);
%%
%Create a posterior predictive samples from the variational model
for i = 1:2
    [mu_r, pcs_r, stddev_r] = variationalBayesSample(imean, psipcs(:, 1:P), psistddev(1:P), lambdapcs(:, 1:P), lambdastddev(1:P), nu);
    %Get a multivariate random sample from the random parameters
    smp(:, i) = multivariateNormalSample(mu_r, pcs_r, stddev_r, 1);    
end
%Display the two samples
hold off
displayShape(smp(:, 1), 'g', 0.7);
hold on
displayShape(smp(:, 2), 'b', 0.7);






