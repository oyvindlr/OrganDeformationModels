[patientData, faces] = createTestPatients;

%%
[popmean, ppcs, pstd, indivmean] = populationModelParameters(patientData);

displayShape(popmean, 'b', 0.4);
hold on;
displayShape(popmean + ppcs(:, 3)*pstd(3), 'g', 0.5);
displayShape(popmean - ppcs(:,3)*pstd(3), 'r', 0.4);

%%
kappa = 0.25;
nu = 6;
[imean, ipcs, istddev] = normalInverseWishartPosteriorParams(vec(patientData(1).contourPoints{2}), popmean, ppcs, pstd, kappa, nu);
displayShape(imean, 'b', 0.4);
hold on;
displayShape(imean + ipcs(:, 1)*istddev(1), 'g', 0.5);
displayShape(imean - ipcs(:, 1)*istddev(1), 'r', 0.4);