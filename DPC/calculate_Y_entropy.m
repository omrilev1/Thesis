clear all; close all; clc;

SNR = -10:1:10;
SNRlin = 10.^(SNR/10);
Px = 1/3;

dx = 1e-3; xAxis = -22:dx:22;

alpha = [0:0.05:1];

uniDist = zeros(size(xAxis)); 
uniDist(abs(xAxis) <=1) = 0.5;

triDist = conv(uniDist,uniDist,'same') * dx;
figure; hold all
h = zeros(length(SNR),length(alpha));
currLegend = [];
for i=1:length(SNR)
    for j=1:length(alpha)
        
        if alpha(j) == 0
            h(i,j) = 1;
        else
            uniAlphaDist = zeros(size(xAxis));
            uniAlphaDist(abs(xAxis) < 1/alpha(j)) = 0.5*alpha(j);
            
            sigma = 1/sqrt(3*SNRlin(i));
            noiseDist = (1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/(sigma^2));
            
            % calculate convolution
            temp = conv(triDist,noiseDist,'same') * dx;
            yDist = conv(temp,uniAlphaDist,'same') * dx;
            
            h(i,j) = log2(alpha(j)) -1*sum(yDist(yDist > 0).*log2(yDist(yDist > 0))) * dx;
        end
    end
    display(strcat('finished SNR = ',num2str(SNR(i))))
    plot(alpha,h(i,:))
end

            
            