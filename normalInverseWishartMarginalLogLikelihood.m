function f = normalInverseWishartMarginalLogLikelihood(data, mu, kappa, Psi_pcs, Psi_stddevs, nu)

J = size(data, 2);
%P = size(data, 1);
P = length(Psi_stddevs);
%P = 1;

if nu <= P+1
    f = -inf;
    return;
end

Psi_stddevs = sqrt(nu-P-1)*Psi_stddevs;

kappastar = kappa + J;
nustar = nu + J;

kappafactor = (P/2)*(log(kappa)-log(kappastar));
gammafactor = multivariateGammaRatioLn(nu/2, nustar/2, P);

m = mean(data, 2);
newcolumns = [data-m sqrt(kappa*J/kappastar)*(m-mu)];
%Project onto the subspace spanned by Psi
newcolumns = Psi_pcs*(Psi_pcs'*newcolumns);

%Dstar = [Psi_pcs*diag(Psi_stddevs) data-m sqrt(kappa*J/kappastar)*(m-mu)];
Dstar = [Psi_pcs*diag(Psi_stddevs) newcolumns];
stddevstar = svd(Dstar, 0);
%stddevstar = stddevstar(1:P);

%A little trick: When computing the pseudo-determinant as the product of the
%non-zero eigenvalues, we want both matrices to have the same number of non-zero
%eigenvalues. We don't want to lose any of the variance of Psi-star, so we
%remove some eigenvalues, but add the corresponding variance to the
%remaining eigenvalues. 
%stddevstar = stddevstar.^2; %Convert to variance/eigenvalues
%stddevstar(1:P) = stddevstar(1:P) + sum(stddevstar(P+1:end))/P;%Add equal parts of removed variance to each remaining eigenvalue
%stddevstar(P) = stddevstar(P) + sum(stddevstar(P+1:end));
%stddevstar = sqrt(stddevstar(1:P));%Remove extra eigenvalues and convert back to standard deviation

psifactor = sum(log(Psi_stddevs))*nu - sum(log(stddevstar(1:P)))*(nu+J);

f = kappafactor + gammafactor + psifactor;




