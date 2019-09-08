function [y] = myModulo(x,delta)
% calculates [x]_{\Delta}. The resuls is in [-Delta,Delta]
y = mod(x + delta,2*delta) - delta;
end

