function H = nnEntropy(samples)
%% nnEntropy calculates the nearest neighbors entropy of samples.
% Input is samples by dimensions

k    = size(samples,2);
n    = size(samples,1);

Ak = (k*pi^(k/2))/gamma(k/2+1);

% Compute nearest neighbor distances
distmat         = squareform(pdist(samples,'euclidean'));
distmat         = distmat + 1e2*eye(n);
md              = min(distmat);

logMD(md > 0)   = log2(md(md > 0));
logMD(md == 0)  = 0;

H = k*mean(logMD) + log2(n*Ak/k) - psi(1)/log(2);
