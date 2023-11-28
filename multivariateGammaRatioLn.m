function g = multivariateGammaRatioLn(a, b, p)
%MULTIVARIATEGAMMARATIOLN Gives the logarithm of the ratio between two 
% multivariate gamma functions with the same dimension p; Gamma_p(b)/Gamma_p(a)
% The difference between the two values, i.e. b-a, must be a multiple of 1/2.

g = 0;

diff = b-a;

for i = 0.5:0.5:diff
    g = g + gammaln(a + i);
end

for i = -(p-1)/2:0.5:(diff-((p-1)/2)-0.5)
    g = g - gammaln(a+i);
end

g2 = 0;

for i = 1:p
    g2 = g2 + gammaln(b + (1 - i)/2);
    g2 = g2 - gammaln(a + (1 - i)/2);
end
g = g2;

% g = log(a);
% for i = 1:diff-1
%     g = g + log(a+i);
% end
% 
% g = p*g;

