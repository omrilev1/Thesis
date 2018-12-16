function [ MI ] = numeric_MI_integral( pdfX,intervalX,pdfY,intervalY,intervalXY,pdfXY )
% This function try to evaluaate numerically the MI between 2 input
% pdfs (X,Y)
% The integral is over (pdfXY*log2(pdfXY/(pdfX*pdfY)))

%% fix intervals , to be in the same resolution 
minX = min(min(intervalX),min(intervalXY,2));
maxX = max(max(intervalX),max(intervalXY,2));
minY = min(intervalY,min(intervalXY,1));
maxY = min(intervalY,max(intervalXY,1));

fixedIntervalX = linspace(minX,maxX,length(pdfX));
fixedIntervalY = linspace(minY,maxY,length(pdfY));
fixedIntervalXY = meshgrid(fixedIntervalX,fixedIntervalY);

dX = mean(diff(fixedIntervalX));
dY = mean(diff(fixedIntervalY));
%% Interpolate pdfs to new intervals
fixedPDFx = interp1(intervalX,pdfX,fixedIntervalX);
fixedPDFy = interp1(intervalY,pdfY,fixedIntervalY);
fixedPDFxy = interp2(intervalXY,pdfXY,fixedIntervalXY);

%% Calculate integrals 
repPDFx = repmat(fixedPDFx,length(fixedPDFy),1);
repPDFy = repmat(fixedPDFy.',1,length(fixedPDFx));

% clculate pdf(x,y)/(pdf(x)*pdf(y))
logTerm = zeros(size(repPDFx));
logTerm(repPDFx > 0 && repPDFy > 0) =  log2(fixedPDFxy./(repPDFx.*repPDFy));

MI = -1*sum(sum(fixedPDFxy.*logTerm)) * dX * dY;

end

