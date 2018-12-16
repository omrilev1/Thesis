%% data
clc;
clear;
snri = (-10:1:24);
snrLin = 10.^(snri/10);
pam = [2,4,8];
capacity = zeros(length(pam),length(snri));

%% Uniform input + PAM input capacity - uniform interferer
InstCapUni = zeros(size(snri));
InstCapDPC = zeros(size(snri));
dX = 0.001;


%% Generate augmented constellation outside the main for loop
for i=1:length(pam)
    Delta = pam(i) + dX;
    xAxis = -8*Delta:dX:8*Delta;
    xAxisForUniformDist = -Delta:dX:Delta;
    
    % inputs - pam
    pamSymbols = (-1*(pam(i) - 1) : 2 : (pam(i) - 1));
    pamOrigDist = ones(size(pamSymbols))/pam(i);
    pamDist = zeros(size(xAxis));
    tempPamIdx = sum(abs(xAxis(:) - pamSymbols) < dX/2,2);
    pamDist(tempPamIdx > 0.5) = 1/pam(i);
    
    % uniform
    uniformDist = zeros(size(xAxis));
    uniformDist(xAxis >= -1*Delta & xAxis <= Delta) = 1/(2*Delta);
    
    % augmented constellation : (v tilde)augmented pam with non uniform probabilities;
    augPamDist = zeros(size(xAxis));
    augPam = [((1 : 2 : (pam(i) - 1)) - 2*Delta) (-1*(pam(i) - 1) : 2 : (pam(i) - 1)) ((-1*(pam(i) - 1) : 2 : -1) + 2*Delta)];
    tempPamIdx = sum(abs(xAxis(:) - augPam) < dX/2,2);
    tempPamProbs = [ abs((-1*(pam(i) - 1) : 2 : -1))/(2*Delta) ...
        (2*Delta - abs((-1*(pam(i) - 1) : 2 : (pam(i) - 1))))/(2*Delta) ...
        abs((-1*(pam(i) - 1) : 2 : -1))/(2*Delta)]/pam(i);
    augPamDist(tempPamIdx > 0.5) = tempPamProbs/dX;
    
    % Calculate composite noise distribution
    oneProbability = zeros(size(pamSymbols));minusOneProbability = zeros(size(pamSymbols));
    
    oneProbability(pamSymbols >= 0) = abs(pamSymbols(pamSymbols >= 0))/(2*Delta);
    minusOneProbability(pamSymbols <= 0) = abs(pamSymbols(pamSymbols <= 0))/(2*Delta);
    probs = [1-abs(pamSymbols(:))/(2*Delta) oneProbability(:) minusOneProbability(:)];
    centers = [0*ones(length(pamSymbols),1) 2*Delta*ones(length(pamSymbols),1) -1*2*Delta*ones(length(pamSymbols),1)];
    
    % calculate original probabilities entropy
    logProbs = log2(probs); logProbs(logProbs == -inf) = 0;
    origConditionedEntropy = sum(logProbs.*probs,2);
    naiveNoiseEntropy = -1*dX*(pamOrigDist * origConditionedEntropy);
    
    
    for j=1:length(snri)
        
        sigma = sqrt((Delta^2)/(3 * 10^(snri(j)/10)));
        %% calculate distributions
        % noise
        noiseDist = 1/sqrt(2*pi*sigma^2) * exp(-0.5*(xAxis.^2)/sigma^2);
        
        % calculate channel output distribution
        yDistFiniteDPC = conv(noiseDist,augPamDist) * dX;
        
        yDist = conv(noiseDist,uniformDist) * dX;
        
        % calculate entropy using gauss-hermite approx. set minimum to original
        % probabilities entropy
        
        noiseConditionedEntropy = gaussianMixtureEntropyApprox(centers,probs,sigma,20);
        
        curNoiseEntropyDPC = -1*sum(noiseConditionedEntropy)/(pam(i));
        
        % numeric integration
        curEntropyFiniteDPC = -1*sum(yDistFiniteDPC(yDistFiniteDPC > 0).*log2(yDistFiniteDPC(yDistFiniteDPC > 0))) * dX;
        curEntropyUniform = -1*sum(yDist(yDist > 0).*log2(yDist(yDist > 0))) * dX;
        
        InstCapUni(i,j) = (curEntropyUniform - 0.5*log2(2*pi*exp(1)*sigma^2));
        
        if isnan(curNoiseEntropyDPC)
            InstCapDPC(i,j) = log2(pam(i));
        else
            InstCapDPC(i,j) = (curEntropyFiniteDPC - curNoiseEntropyDPC);
        end
    end
end
% AWGN Channel capacity
GaussianCapacity = 0.5*log2(1+ 10.^(snri/10));


%% Plot Graphs - as function of SNR and Eb/N0
figure;hold all
plot(snri,GaussianCapacity,'LineWidth',1.6);
plot(snri,InstCapUni,'LineWidth',1.6);
plot(snri,InstCapDPC,'LineWidth',1.6);
set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('SNR (Es/No) dB','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');
legend('Gaussian Capacity','Gaussian Channel Uniform Input (Regular DPC)',...
    'PAM2 , uniform interferer I(v;Y)',...
    'PAM4 , uniform interferer I(v;Y)',...
    'PAM8 , uniform interferer I(v;Y)');
title('Capacity for AWGN , With Uniform Input and uniform interferer')

figure;hold all
plot(10*log10(snrLin./(2*GaussianCapacity)),GaussianCapacity,'LineWidth',1.6);
plot(10*log10(snrLin./(2*InstCapUni)),InstCapUni,'LineWidth',1.6);
plot(10*log10(snrLin(2:end)./(2*InstCapDPC(2:end))),InstCapDPC(2:end),'LineWidth',1.6);
set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('(Eb/No) dB','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');
legend('Gaussian Capacity','Gaussian Channel Uniform Input (Regular DPC)',...
    'Finite Power DPC , PAM input I(v;Y)')
title('Capacity for AWGN , With Uniform Input and uniform interferer')
