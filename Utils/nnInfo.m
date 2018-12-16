function I = nnInfo(x,y)
%% nnInfo calculates mutual information between
% x and y by using nearest neighbors entropy and
% I = HX + HY - HXY

I = nnEntropy(x) + nnEntropy(y) - nnEntropy([x,y]);
