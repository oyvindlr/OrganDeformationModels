%Create some data to run examples with
[patientData, faces] = createTestPatients;

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
%The parameters for the population model are also used in the NIW-model.
%Here, we need additional parameters kappa and nu
kappa = 0.25;
nu = 6;
% Calculate the posterior parameters given some patient data.
% The posterior parameters can also be used as point estimates for the mean
% and covariance matrix (represented here by principal components) for that patient
[imean, ipcs, istddev] = normalInverseWishartPosteriorParams(vec(patientData(1).contourPoints{2}), popmean, ppcs, pstd, kappa, nu);

%Show the mean shape and +/-one standard deviation of the first mode
displayShape(imean, 'b', 0.4);
hold on;
displayShape(imean + ipcs(:, 1)*istddev(1), 'g', 0.5);
displayShape(imean - ipcs(:, 1)*istddev(1), 'r', 0.4);
hold off
%%
%Calculate the parameters for a variational Bayes model
[popmean, pcs_psi, stddev_psi, pcs_lambda, stddev_lambda] = variationalBayesModelParameters(patientData, nu);

%Show the mean shape and +/-one standard deviation of the first INTER-PATIENT mode
displayShape(popmean, 'b', 0.4);
hold on;
displayShape(popmean + pcs_lambda(:, 1)*stddev_lambda(1), 'g', 0.5);
displayShape(popmean - pcs_lambda(:,1)*stddev_lambda(1), 'r', 0.4);
hold off

%%
%Calculate the patient-specific posterior parameters of the variational
%Bayes model. We need to assign a few extra parameter values
delta_psi = 1;
delta_lambda = 1;
D_psi = pcs_psi*diag(stddev_psi);
D_lambda = pcs_lambda*diag(stddev_lambda);
patData = [vec(patientData(1).contourPoints{1}) vec(patientData(1).contourPoints{2})];

[imean, psipcs, psistddev, psireg, lambdapcs, lambdastddeev, lambdareg] = ...
    fastVariational(popmean, D_psi, D_lambda, delta_psi, delta_lambda, patData, nu);

%Display the expected mean shape
displayShape(imean);

%%
%Create a new patient from the variational model
[mu, pcs_n, stddev_n] = newPatientDistributionVariational(popmean, pcs_psi, stddev_psi, pcs_lambda, stddev_lambda, nu);
