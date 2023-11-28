function c = niw_costfunc(avector, patientData, popmean, ipcs, istddev, J)
     kappa = avector(1);
     nu = avector(2);
     penalty = 0;
     p = length(istddev);
     %force nu to be greater than or equal to P, otherwise the function
     %cannot be evaluated
     if nu <= p + 1 + 1e-2
         %add a penalty, so the optimization needs to incrase nu
         penalty = 10*(p + 1 + 1e-2 - nu)^2;
         
         nu = p + 1 + 1e-2;         
     end
     if kappa <= + 1e-2
         penalty = penalty + 10*(kappa - 1e-2)^2;
         kappa = 1e-2;
     end
     
     c = -totalLogLikelihood(patientData, popmean, kappa, ipcs, istddev, nu, J);
     c = c + penalty;
end