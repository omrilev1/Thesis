function H = nnEntropy(samples)
%% nnEntropy calculates the nearest neighbors entropy of samples.
% Input is samples by dimensions

k    = 1;
n    = length(samples);

Ak = (k*pi^(k/2))/gamma(k/2+1);

% Compute nearest neighbor distances
distmat         = squareform(pdist(samples(:),'euclidean'));
distmat         = distmat + 1e2 * eye(n);
md              = max(distmat);

logMD(md > 0)   = log2(md(md > 0));
logMD(md == 0)  = 0;

H = k*mean(logMD) + log2(n*Ak/k) - psi(1)/log(2);
