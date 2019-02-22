function [f,cordX,cordY,dx] = pdfGenerator(Type,variance,n)
% This Function gets as input pdf parameters (type , variance and
% dimension) and output the appropriate pdf

if strcmp(Type, 'Gaussian')
    
    baseCord = linspace(-8*sqrt(variance),8*sqrt(variance),8e2);
    dx = 2*8*sqrt(variance) / 8e2;
    [gridX,gridY] = meshgrid(baseCord,baseCord);
    
    tempPDF = dx*(1/sqrt(2*pi*variance))^2 * exp(-1*(gridX.^2 + gridY.^2)/(2*variance));
    
elseif  strcmp(Type, 'Exp')
    
    baseCord = linspace(0,8*sqrt(variance),8e2);
    dx = 8*sqrt(variance) / 8e2;
    [gridX,gridY] = meshgrid(baseCord,baseCord);
    
    lambda = 1/sqrt(variance);
    tempPDF = dx*(lambda)^2 * exp(-1*(lambda*gridX + lambda*gridY));
    
elseif  strcmp(Type, 'Laplace')
    
    baseCord = linspace(-8*sqrt(variance),8*sqrt(variance),8e2);
    dx = 2*8*sqrt(variance) / 8e2;
    [gridX,gridY] = meshgrid(baseCord,baseCord);
    
    b = sqrt(variance/2);
    tempPDF = dx*(((1/(2*b)))^2)* exp(-1*(abs(gridX)/b + abs(gridY)/b));    
elseif strcmp(Type, 'Unifrom')
    
    baseCord = linspace(-8*sqrt(variance),8*sqrt(variance),8e2);
    dx = 2*8*sqrt(variance) / 8e2;
    [gridX,gridY] = meshgrid(baseCord,baseCord);
    
    Delta = sqrt(3*variance);
    tempPDF = zeros(size(gridX)); tempPDF((gridX > -Delta)&&(gridX < Delta) && (gridY > -Delta)&&(gridY < Delta)) = dx*(1/(2*Delta))^2;
    
end

if n==1
    f = sum(tempPDF,1);
    cordX = gridX(1,:);
    cordY = [];
else
    f = tempPDF;
    cordX = gridX;
    cordY = gridY;
end

end

