function [mu0, psi_pcs, psi_stddevs, psi_regularization, lambda_pcs, lambda_stddevs, lambda_regularization] = fastVariational(mu0, DL, DP, dL, dP, S, nu)

%FASTVARIATIONAL Implements the computationally efficient computation of
%posterior parameters for the variational Bayes organ deformation model 
%described in Rørtveit et al., 2023
%
% function [mu0, psi_pcs, psi_stddevs, psi_regularization, lambda_pcs, 
%                            lambda_stddevs, lambda_regularization]
%                             = fastVariational(mu0, DL, DP, dL, dP, S, nu)
%
%
% Input arguments (see paper for details):
%  mu0: Mean vector
%  DL:  Data matrix for Lambda - the inter-patient covariance matrix
%  DP:  Data matrix for Psi - a multiple of the intra-patient covariance
%       matrix
%  dL:  delta_Lambda - regularization coefficient for Lambda
%  dP:  delta_Psi - regularization  coefficient for Psi
%  S:   Matrix containing the patient specific data. Each column represents
%       one observed shape vector
%  nu:  Degrees of freedom parameter of the inverse-Wishart
%       distribution
%
% Output arguments:
%  mu0:             Posterior mean vector
%  psi_pcs:         Principal components of the posterior Psi-matrix
%  psi_stddevs:     Standard deviations of the posterior Psi-matrix PCs
%  psi_regularization: Regularization factor of the posterior Lambda-matrix
%
%                      The posterior Psi is equal to
%                      Psi_pcs*psi_stddevs*psi_stddevs*Psi_pcs'+psi_regularization*eye(P)
%
%  lambda_pcs:         Principal components of the posterior Lambda-matrix
%  lambda_stddevs:     Standard deviations of the posterior Lambda-matrix PCs
%  lambda_regularization: Regularization factor of the posterior Lambda-matrix


[P, n] = size(S);
NL = size(DL, 2);
NP = size(DP, 2);

DP = sqrt(nu-NP-1)*DP;
dP = (nu-NP-1)*dP;


%%% The commented code here and in the loop computes the values of the
%%% various parameters as they would be if we did not use the
%%% fast-computation tricks. This way, the fast computation can be compared
%%% to the much simpler and slower direct computation of the same values.
%%%%%
% Psi = DP*DP'+dP*eye(P);
% PsiI = Psi;
% Lambda = DL*DL'+dL*eye(P);
% mu0_f = mu0;
% qk = Lambda\mu0;
%%%%%

eyeNP = eye(NP+n);
K = blkdiag(zeros(NL), eyeNP);
Q = blkdiag(-(1/dL)*inv(dL*eye(NL)+(DL'*DL)), zeros(NP+n));
q = (1/dL)*(mu0 - DL*((dL*eye(NL)+(DL'*DL))\(DL'*mu0)));
DPorig = DP;
DP = [DPorig zeros(P, n)];
D = [DL DP];
L = blkdiag(zeros(NL), -(1/dP)*inv(dP*eyeNP+DP'*DP));
dPorig = dP;
lastmu = inf;
i = 0;
s = mean(S, 2);
%nu  = nu - NP -1;
while norm(mu0-lastmu) > 1e-4 && i < 100000
    lastmu = mu0;
    i = i + 1;    
    d = 1/((1/dL)+n*(nu+n)*(1/dP));
    F = Q + n*(n+nu)*L;
    H = inv(D'*D+inv(d*F));
    r = q + n*(n+nu)*((1/dP)*s+D*(L*(D'*s)));
    mu0 = d*r - d*D*(H*(D'*r));
    %%%%%%
%     LambdaIfast = d*eye(P)-d*D*H*D';
%     LambdaI = inv(inv(Lambda)+n*(nu+n)*inv(PsiI));
%     ra = qk + n*(nu+n)*(PsiI\s);
%     mu0_f = LambdaI*ra;
    %%%%%%
    DP = [DPorig S-mu0];
    D = [DL DP];
    dP = dPorig + n*d;
    G = K-n*d*H;
    L = -(1/dP)*inv(dP*inv(G)+D'*D);
    %%%%%%
%     PsiI = Psi + (S-mu0_f)*(S-mu0_f)' + n*LambdaI;
%     PsiIfast = dP*eye(P)+D*G*D';
    %%%%%
end
%nu = nu + NP + 1;

%PCA
[psi_pcs, psi_stddevs] = svd(D*(chol(G)'), 0);
psi_stddevs = diag(psi_stddevs);
psi_stddevs = 1/sqrt(nu-NP-1)*psi_stddevs;
psi_regularization = 1/(nu-NP-1)*dP;
[lambda_pcs, lambda_stddevs] = svd(sqrt(d)*D*(chol(-H)'), 0);
lambda_stddevs = diag(lambda_stddevs);
lambda_regularization = d;




