function [ outputPDF,outputInterval ] = pdfAfterConstAddition( inputPDF,inputInterval,mu )
% This function calculate the pdf of variable + mu 
% (mu can be vector and then we'll have 2 dimensional distribution). 
% for example , if we have Y = X + mu , then the input pdf will be the pdf
% of X , and the output will be the pdf of Y
% The mathematical relation is pdf(x-mu)

outputInterval = inputInterval(:) + (mu(:)).';
outputPDF = inputPDF;

end