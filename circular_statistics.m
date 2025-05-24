function [R, Z, P, M, K] = circular_statistics(phiv, weights)
N = length(phiv);
if nargin < 2
    sinm = mean(sin(phiv));
    cosm = mean(cos(phiv));
else
    sinm = sum(sin(phiv).*(weights/sum(weights, 'omitnan')), 'omitnan');
    cosm = sum(cos(phiv).*(weights/sum(weights, 'omitnan')), 'omitnan');
end
R = sqrt(sinm^2+cosm^2);
Z = N * R^2;
if N <= 50
    P = exp(-Z)*(1+(2*Z-Z^2)/(4*N)-(24*Z-132*Z^2+76*Z^3-9*Z^4)/(288*N^2));
else
    P = exp(-Z);
end
M = atan2(sinm, cosm);
K = R*(2-R^2)/(1-R^2);
end