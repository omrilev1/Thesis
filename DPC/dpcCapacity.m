%% Evaluate the MI in various DPC problems
close all; clear; clc;

addpath(genpath('Utils'));
% regular DPC : compare between the regular AWGN setting (awgn with
% gaussian input) , AWGN with uniform input and AWGN with interferer and
% DPC
calcMethod = 'nnEntropy'; %'semiAnalytic' for working directly on the pdfs
% define parameters
Delta = 2; % lattice constant , 1D
Px = Delta^2 / 12;

SNR = -4:2:24; %-10:2:45;
SNRlin = 10.^(SNR/10);
sigmaZ = sqrt(Px./SNRlin);

NforStatistics = 1e4;
NforMI = 1e3; iter = 1;
IxyAWGN        = zeros(size(SNR));
IxyAWGNuni = zeros(size(SNR));
IxyDPC         = zeros(size(SNR));

for i=1:length(SNR)
    
    % generate random stream  : uniform between [0,Delta] and wgn
    awgnInput = sqrt(Px) * randn(1,NforStatistics);
    uniInput  = Delta * rand(1,NforStatistics);
    
    % generate noise : wgn (for regular awgn and uniform input chanels) and
    % alpha*wgn + (1-alpha)*uniform modulo for DPC
    Z = sigmaZ(i)*randn(1,NforStatistics);
    U = Delta * rand(1,NforStatistics);
    
    
    % generate modulo noise with optimal scaling factor
    alpha = SNRlin(i)/(SNRlin(i) + 1);
    Zprime = (1-alpha)*U + alpha*Z;
    
    ZprimeMod = mod(Zprime,Delta);
    
    % generate channel output
    yAWGN = awgnInput + Z;
    yAWGNuniform = uniInput + Z;
    yDPC = uniInput + ZprimeMod;
    yDPC = mod(yDPC,Delta);
    
    % Calculate entropies
    switch calcMethod
        case 'semiAnalytic'
            % semi analytic
            dX = 0.01;
            x = -50:dX:50;
            GaussianInputPDF = (1/sqrt(2*pi*Px))*exp(-x.^2/(2*Px));
            GaussianNoisePDF = (1/sqrt(2*pi*sigmaZ^2))*exp(-x.^2/(2*sigmaZ^2));
            uniInputPDF = zeros(size(x)); uniInputPDF(x >=-Delta/2 && x<= Delta/2) = 1/Delta;
            
            % calculate modulo noise pdf 
            scaledNoisePDF = (1/sqrt(2*pi*alpha^2*sigmaZ^2))*exp(-x.^2/(2*alpha^2*sigmaZ^2));
            scaledUniPDF = zeros(size(x)); scaledUniPDF((1-alpha)*x >=0 && x<= Delta/(1-alpha)) = (1-alpha)/Delta;
            
            zPrimePDF = conv(scaledNoisePDF,scaledUniPDF);
            zPrimeSupport = linspace(-2*min(x),2*max(x),length(zPrimePDF));
            
            zPrimeModuloSupport = linspace(0,Delta,length(zPrimeSupport));
            zPrimeModuloPDF = zeros(size(zPrimeModuloSupport));
            
            
        otherwise
            % numeric
            [HxAWGN,HyAWGN,HxyAWGN,IxyAWGN(i)] = mutualInfoCalculator(awgnInput,yAWGN,NforMI,iter,calcMethod);
            [HxAWGNuni,HyAWGNuni,HxyAWGNuni,IxyAWGNuni(i)] = mutualInfoCalculator(uniInput,yAWGNuniform,NforMI,iter,calcMethod);
            [HxDPC,HyDPC,HxyDPC,IxyDPC(i)] = mutualInfoCalculator(uniInput,yDPC,NforMI,iter,calcMethod);
    end
    
    
    disp(strcat('Finished iteration for SNR = ',num2str(SNR(i))))
end
awgnCapTheory = 0.5*log2(1 + SNRlin);
figure;
plot(SNR,awgnCapTheory,'--','LineWidth',1.5);hold on;
plot(SNR,IxyAWGN,'-*','LineWidth',1.5);hold on;
plot(SNR,IxyAWGNuni,'-o','LineWidth',1.5);hold on;
plot(SNR,IxyDPC,'-^','LineWidth',1.5);hold on;
grid on; grid minor;
xlabel('SNR [dB]'); ylabel('I(X;Y)');
title('Capacity in various scenarios')
legend('Theoretic AWGN','Simulation AWGN','AWGN with Uniform input','DPC');

% Capacity(Eb/N0)
figure;
plot(10*log10(SNRlin./(2*awgnCapTheory)),awgnCapTheory,'--','LineWidth',1.5);hold on;
plot(10*log10(SNRlin./(2*IxyAWGN)),IxyAWGN,'-*','LineWidth',1.5);hold on;
plot(10*log10(SNRlin./(2*IxyAWGNuni)),IxyAWGNuni,'-o','LineWidth',1.5);hold on;
plot(10*log10(SNRlin./(2*IxyDPC)),IxyDPC,'-^','LineWidth',1.5);hold on;
grid on; grid minor;
xlabel('E_b/N_0 [dB]'); ylabel('I(X;Y)');
title('Capacity in various scenarios')
legend('Theoretic AWGN','Simulated AWGN','AWGN With Uniform Input','DPC With Modulo');
