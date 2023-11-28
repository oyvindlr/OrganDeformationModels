function l = multivariateNormalLogLikelihood(x,mu, pcs, stddevs)


stddevinv = 1./stddevs;

rfac = -sum(log(stddevs));

d = diag(stddevinv)*pcs'*(x-mu);
d = d.^2;

l = size(x, 2)*rfac - 1/2*sum(sum(d));
