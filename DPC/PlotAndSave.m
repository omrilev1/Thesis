% This script calculate capacity for different DPC settings. The main 
% purpose is to calculate the capacity with finite power interferer ,
% without modulo in the reciever.
clc; clear ; close all;
%% data
snri = (-10:2:20);
snrLin = 10.^(snri/10);
pam = 2; % [2,4,8,16,32,64];
dX = 0.001; % for numeric calculations

%% Uniform input + PAM input capacity - uniform interferer 
InstCapUni = zeros(size(snri));
InstCapDPCbound = zeros(size(snri));
InstCapDPC = zeros(size(snri));

capacityPAM = zeros(length(pam),length(snri));
capacityUNI = zeros(1,length(snri));

for i=1:length(snri)
    Delta = sqrt(3*10^(snri(i)/10));
    xAxis = -5*Delta:dX:5*Delta;
    xAxisForUniformDist = -Delta:dX:Delta;
    
    %% calculate distributions
    % noise
    noiseDist = 1/sqrt(2*pi) * exp(-0.5*(xAxis.^2));
    
    % inputs - uniform and pam
    uniformDist = zeros(size(xAxis));
    uniformDist(xAxis >= -1*Delta & xAxis <= Delta) = 1/(2*Delta);
    
    % augmented constellation
    triangleDist = conv(uniformDist,uniformDist)*dX;
    ValidIdx = round(length(uniformDist)/2):(length(triangleDist)-round(length(uniformDist)/2));
    triangleDist = triangleDist(ValidIdx);
        
    % calculate channel output distribution
    yDist = conv(noiseDist,uniformDist) * dX;
    yDistFiniteDPC = conv(noiseDist,triangleDist) * dX;
    
    % Calculate composite noise distribution
    oneProbability = zeros(size(xAxisForUniformDist));minusOneProbability = zeros(size(xAxisForUniformDist));
    
    oneProbability(xAxisForUniformDist >= 0) = abs(xAxisForUniformDist(xAxisForUniformDist >= 0))/(2*Delta);
    minusOneProbability(xAxisForUniformDist <= 0) = abs(xAxisForUniformDist(xAxisForUniformDist <= 0))/(2*Delta);
    probs = [1-abs(xAxisForUniformDist(:))/(2*Delta) oneProbability(:) minusOneProbability(:)];
    centers = [0*ones(length(xAxisForUniformDist),1) 2*Delta*ones(length(xAxisForUniformDist),1) -1*2*Delta*ones(length(xAxisForUniformDist),1)];
    
    % calculate entropy using gauss-hermite approx. set minimum to original
    % probabilities entropy
    logProbs = log2(probs); logProbs(logProbs == -inf) = 0;
    origConditionedEntropy = sum(logProbs.*probs,2);
    naiveNoiseEntropy = -1*dX*sum(origConditionedEntropy)/(2*Delta) + 0.5*log2(2*pi*exp(1));
    
    noiseConditionedEntropy = gaussianMixtureEntropyApprox(centers,probs,1,20);
    
    curNoiseEntropyDPC = -1*dX*sum(noiseConditionedEntropy)/(2*Delta);
    
    % numeric integration
    curEntropy = -1*sum(yDist(yDist > 0).*log2(yDist(yDist > 0))) * dX;
    curEntropyFiniteDPC = -1*sum(yDistFiniteDPC(yDistFiniteDPC > 0).*log2(yDistFiniteDPC(yDistFiniteDPC > 0))) * dX;
    
    InstCapUni(i) = (curEntropy - 0.5*log2(2*pi*exp(1)));
    InstCapDPCbound(i) = (curEntropyFiniteDPC - 0.5*log2(2*pi*exp(1)));
    
    if isnan(curNoiseEntropyDPC)
        InstCapDPC(i) = (curEntropyFiniteDPC - naiveNoiseEntropy);
    else
        InstCapDPC(i) = (curEntropyFiniteDPC - curNoiseEntropyDPC);
    end
end
% AWGN Channel capacity
GaussianCapacity = 0.5*log2(1+ 10.^(snri/10));
% Finite Power DPC Upper Bound
FiniteDPCbound = 0.5*log2(1+ 2*10.^(snri/10));

%% Plot Graphs - as function of SNR and Eb/N0
figure;hold all
plot(snri,GaussianCapacity,'LineWidth',1.6);
plot(snri,FiniteDPCbound,'LineWidth',1.6);
plot(snri,InstCapUni,'LineWidth',1.6);
plot(snri,InstCapDPCbound,'LineWidth',1.6);
plot(snri,InstCapDPC,'LineWidth',1.6);
set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('SNR (Es/No) dB','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');
legend('Gaussian Capacity','Finite Power DPC Bound','Gaussian Channel Uniform Input (Regular DPC)',...
    'Finite Power DPC uniform interferer and input I(v\tilde ; Y)','Finite Power DPC I(v;Y)')
title('Capacity for AWGN , With Uniform Input and uniform interferer')

figure;hold all
plot(10*log10(snrLin./(2*GaussianCapacity)),GaussianCapacity,'LineWidth',1.6);
plot(10*log10(snrLin./(2*FiniteDPCbound)),FiniteDPCbound,'LineWidth',1.6);
plot(10*log10(snrLin./(2*InstCapUni)),InstCapUni,'LineWidth',1.6);
plot(10*log10(snrLin./(2*InstCapDPCbound)),InstCapDPCbound,'LineWidth',1.6);
plot(10*log10(snrLin(2:end)./(2*InstCapDPC(2:end))),InstCapDPC(2:end),'LineWidth',1.6);
set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('(Eb/No) dB','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');
legend('Gaussian Capacity','Finite Power DPC Bound','Gaussian Channel Uniform Input (Regular DPC)',...
    'Finite Power DPC uniform interferer and input I(v\tilde ; Y)','Finite Power DPC I(v;Y)')
title('Capacity for AWGN , With Uniform Input and uniform interferer')
