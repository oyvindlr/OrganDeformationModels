function c = niw_point_costfunc(avector, patientData, mu0, pcs, stddevs, N)
     kappa = avector(1);
     penalty = 0;
     nu = avector(2);
     p = length(stddevs);
     if nu <= p + 1 + 1e-2
         %add a penalty, so the optimization needs to incrase nu
         penalty = 100*(p + 1 + 1e-2 - nu)^2;
         
         nu = p + 1 + 1e-2;         
     end
     if kappa <=  1e-2
         penalty = penalty + 100*(kappa- 1e-2)^2;
         kappa = 1e-2;
     end
     c = -niwPosteriorPointestimatesLikelihood(patientData, mu0, kappa, pcs, stddevs, nu, N);
     c = c + penalty;
end