[patientData, faces] = createTestPatients;

%%
[popmean, ppcs, pstd, indivmean] = populationModelParameters(patientData);

displayShape(popmean, 'b', 0.4);
hold on;
displayShape(popmean + ppcs(:, 1)*pstd(1), 'g', 0.5);
displayShape(popmean - ppcs(:,1)*pstd(1), 'r', 0.4);

%%
kappa = 0.25;
nu = 6;
[imean, ipcs, istddev] = normalInverseWishartPosteriorParams(vec(patientData(1).contourPoints{2}), popmean, ppcs, pstd, kappa, nu);
displayShape(imean, 'b', 0.4);
hold on;
displayShape(imean + ipcs(:, 2)*istddev(2), 'g', 0.5);
displayShape(imean - ipcs(:, 2)*istddev(2), 'r', 0.4);
%%
[popmean, pcs_intra, stddev_intra, pcs_inter, stddev_inter] = variationalBayesModelParameters(patientData);
displayShape(popmean, 'b', 0.4);
hold on;
displayShape(popmean + pcs_inter(:, 1)*stddev_inter(1), 'g', 0.5);
displayShape(popmean - pcs_inter(:,1)*stddev_inter(1), 'r', 0.4);


%%
delta_psi = 240000;
delta_lambda = 80000;
alpha = 0.4;
[imean, psipcs, psistddev, psireg, lambdapcs, lambdastddeev, lambdareg] = ...
    fastVariational(popmean, sqrt(nu)*pcs_intra*diag(stddev_intra), ...
            sqrt(alpha)*pcs_inter*diag(stddev_inter), delta_psi, delta_lambda, ...
                    vec(patientData(1).contourPoints{2}), nu);
                
                