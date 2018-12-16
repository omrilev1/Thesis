function [ pdf,interval ] = pdfOfSum( pdf1,interval1,pdf2,interval2 )
% This function gets 2 pdfs of 2 random variables and return the pdf of the
% sum. The function actually calculate the convolution between the
% variables
% pdf1/pdf2           : vector of pdf values
% interval1/interval2 : values pdf evaluated in. Should be a vector of
% the form xStart:dX:xEnd
% The function interpolate the pdfs to be in the same resolution , and then
% calculate the convolution and the resulting output interval

%% calculate the dX of the resulting pdf - will be the maximum , and then
%% we'll calculate the final pdf according to this interval
delta1 = mean(diff(interval1));
delta2 = mean(diff(interval2));

%% Fix pdfs to be with the same resolution
if (delta1 == delta2)
    fixPdf1 = pdf1;
    fixPdf2 = pdf2;
    fixInterval1 = interval1;
    fixInterval2 = interval2;
elseif (delta1 < delta2)
    fixPdf2 = pdf2;
    fixInterval2 = interval2;
    fixInterval1 = min(interval1):delta2:max(interval1);
    fixPdf1 = interp1(interval1,pdf1,fixInterval1);
else
    fixPdf1 = pdf1;
    fixInterval1 = interval1;
    fixInterval2 = min(interval2):delta1:max(interval2);
    fixPdf2 = interp1(interval2,pdf2,fixInterval2);
end

pdf = conv(fixPdf1,fixPdf2);
interval = linspace(min(fixInterval1) + min(fixInterval2),max(fixInterval1) + max(fixInterval2),length(pdf));


end

