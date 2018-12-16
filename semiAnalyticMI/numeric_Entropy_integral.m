function [Entropy] = numeric_Entropy_integral( pdf,interval)
% This function try to evaluate numerically the Entropy of pdf=
% The integral is over (pdfX*log2(pdfX))
dY = mean(diff(interval,1,1),1);

logTerm = zeros(size(pdf));
logTerm(pdf > 0) =  log2(pdf(pdf > 0));

Entropy = -1*sum(pdf.*logTerm,1) * dY;

end

