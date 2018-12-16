function [ outputPDF,outputInterval ] = pdfAfterScaling( inputPDF,inputInterval,alpha )
% This function calculate the pdf of scaled variable. 
% for example , if we have Y = alpha*X , then the input pdf will be the pdf
% of X , and the output will be the pdf of Y
% The mathematical relation is 
% (1/alpha)*pdf(x/alpha)

outputInterval = inputInterval/alpha;
outputPDF = inputPDF/alpha;

end

