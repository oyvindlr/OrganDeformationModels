function l = niwPosteriorPointestimatesLikelihood(patientData, mu0, kappa, pcs, stddevs, nu, N)

l =0;
%neigs = min(N, length(stddevs));
neigs = length(stddevs);
for i = 1:length(patientData)
    d = dmat(patientData(i).contourPoints);
    %m = mean(d(:, 1:N), 2);
    %d = pcs*(pcs'*(d-mu0))+mu0;
    [patmean, pcs_pat, stddev_pat] = normalInverseWishartPosteriorParams(d(:, 1:N), mu0, pcs, stddevs, kappa, nu);
    %pcs_pat = pcs_pat(:, 1:neigs+N);
    %stddev_pat = stddev_pat(1:neigs+N);
    l = l + multivariateNormalLogLikelihood(d(:, N+1:end), patmean, pcs_pat, stddev_pat);     
end

end


function d = dmat(points)
d = zeros(numel(points{1}), length(points));
for i = 1:length(points)
    d(:, i) = vec(points{i});
end
end