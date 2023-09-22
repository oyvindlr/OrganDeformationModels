function D_sample = wishartSample(D, reg, nu)

D_sample = D*randn(size(D, 2), nu) + sqrt(reg)*randn(size(D,1), nu); 