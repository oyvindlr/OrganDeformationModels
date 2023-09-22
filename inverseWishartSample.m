function [pcs_sample, stddev_sample] = inverseWishartSample(psi_pcs, psi_stddev, nu)

invD = psi_pcs*diag(1./psi_stddev);

sample = invD*randn(length(psi_stddev), nu);

[pcs_sample, stddev_inv] = svd(sample, 0);


stddev_sample = 1./diag(stddev_inv);

stddev_sample(size(psi_pcs, 2)+1:end) = 0;

[stddev_sample, ind] = sort(stddev_sample, 'desc');
pcs_sample = pcs_sample(:, ind);


% function D_sample = inverseWishartSample(D, reg, nu)
% 
% M = size(D, 2);
% %Psi inverse = 1/regI - 1/reg*D*(1/reg I+D^TD)^-1 DT
% reg_inv = 1/reg;
% %DItta e gale  minus gange minus blir pluss
% D_inv = -sqrt(reg_inv)*D*(chol(inv(reg*eye(M)+D'*D))');
% 
% D_wish = wishartSample(D_inv, reg_inv, nu);
% 
% 
% [pcs, stds] = svd(D_wish, 0);
% 
% D_sample = pcs*diag(1./diag(stds));
