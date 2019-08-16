function [Hout] = parityMatrix_toSystematic(Hin)
% this function convert parity check matrix to systematic form, i.e
% Hout = [A I_(n-k)]
N = size(Hin,2);
K = -1*(size(Hin,1)-N);

HH = mod(rref(Hin),2);
Hout = [HH(:,(N-K+1):end) HH(:,1:(N - K))];

end

