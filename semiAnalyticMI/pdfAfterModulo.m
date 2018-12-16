function [ moduloPdf,moduloInterval ] = pdfAfterModulo( pdf,interval,moduloConst )
% This function gets pdf of random variable and return the pdf after calculating 
% modulo (modulo interval is also an input). 
% The function "wrappes" the pdf to the output interval , and calculate the
% right interval resolution (well end with the same length pdf, thus the
% resolution will be finer after the modulo)

% pdf   : vector of pdf values
% interval : values pdf evaluated in. Should be a vector of
% the form xStart:dX:xEnd
% moduloConst : we perform modulo to the interval [-moduloConst/2 , moduloConst/2]

% enhance input pdf , to include integer multiple of the basic modulo
% interval
leftEdge = moduloConst * floor(min(interval)/moduloConst) - moduloConst/2;
rightEdge = moduloConst * ceil(max(interval)/moduloConst) + moduloConst/2;

% enhanced pdf has +10% points relative to input pdf (arbitrary choice)
enhancedInterval = linspace(leftEdge,rightEdge,1.1*length(interval));
enhancedPDF = zeros(size(enhancedInterval));
enhancedPDF(enhancedInterval >= min(interval) && enhancedInterval <= max(interval)) = ...
    interp1(interval,pdf,enhanncedInterval(enhancedInterval >= min(interval) && enhancedInterval <= max(interval)));

% Wrapp the pdf to the basic interval
moduloPdf = zeros(1,length(enhancedInterval));
moduloInterval = linspace(-moduloConst/2,moduloConst/2,length(moduloPdf));

for i= (leftEdge + moduloConst/2)/moduloConst : (rightEdge - moduloConst/2)/moduloConst
    curIntervalValues = linspace(i*moduloConst - moduloConst/2,i*moduloConst + moduloConst/2,...
        length(moduloPdf));
    
    % find relevant indices of enhanced PDF
    indices = ((enhancedInterval >= min(curIntervalValues)) &&  (enhancedInterval <= max(curIntervalValues)));
    
    % sum the appropriate ppdf values (perform interpolation if needed)
    moduloPdf = moduloPdf + interp1(enhancedInterval(indices),enhancedPDF(indices),curIntervalValues); 
end

