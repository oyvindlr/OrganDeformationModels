function l = totalLogLikelihood(patientData, mu0, kappa, Psi_pcs, Psi_stddevs, nu, J)


l= 0;
P = numel(patientData(1).contourPoints{1});
for i = 1:length(patientData)
    K = min(J, length(patientData(i).contourPoints));
    dataMat = zeros(P, K);
    for j = 1:K
        dataMat(:, j) = vec(patientData(i).contourPoints{j});
    end
    %m = mean(dataMat, 2);
    %dataMat = Psi_pcs*(Psi_pcs'*(dataMat - mu0)) + mu0;
    l = l + normalInverseWishartMarginalLogLikelihood(dataMat, mu0, kappa, Psi_pcs, Psi_stddevs, nu);
end
        
        