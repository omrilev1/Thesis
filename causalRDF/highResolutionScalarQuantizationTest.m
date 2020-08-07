
clear all; clc; close all

% uniform quantizer
delta = 0.025;% 1.5/sqrt(2);
type = 'centroid';
partition = (-6.5:delta:6.5);
switch type
    case 'regular'
        codebook = (-7:1:7)*delta;
    case 'centroid'
        [codebook,~] = generateGaussianCodeBook(partition);
        %probs = [logspace(log(1e-3),log(1.15),(length(codebook) - 1)/2) 1 fliplr(logspace(log(1e-3),log(1.15),(length(codebook) - 1)/2))];
        probs = [logspace(log(1e-3),log(1.15),(length(codebook) - 1)/2) 1 fliplr(logspace(log(1e-3),log(1.15),(length(codebook) - 1)/2))];
end
centers = [partition(1) (partition(1:end-1) + partition(2:end))/2 partition(end)];

x = randn(1,4*1e4);
currDistance = abs(x - centers(:));
[~,minIdx] = min(currDistance,[],1);
T = codebook(minIdx);
S = x - T;
figure; histogram(S);
title(strcat('Quantization Error: mean = ',num2str(mean(S)),'Var/(\Delta^2/12) = ',num2str((mean(S.^2) - (mean(S))^2)/(delta^2/12))))
grid on; grid minor;

function [codebook,probs] = generateGaussianCodeBook(partition)

dx = 0.00001;
x = (-3 + partition(1)):dx:(partition(end) + 3);

pdf = (1/sqrt(2*pi))*exp(-0.5*x.^2);
codebook = zeros(1,length(partition) - 1);
probs = zeros(1,length(partition) - 1);

for i=1:(length(partition) - 1)
    
    curIdx = find((x >= partition(i)) & (x <= partition(i+1)));
    curPDF = pdf(curIdx) / sum(pdf(curIdx));
    curX = x(curIdx);
    codebook(i) = sum(curX.*curPDF);
    probs(i) = sum(pdf(curIdx))/sum(pdf);
    
end
end