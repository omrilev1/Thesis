addpath(genpath('..\Utils'))

% This script calculate capacity for different DPC settings. The main
% purpose is to calculate the capacity with finite power interferer ,
% without modulo in the reciever.
clc; clear ; close all;
%% data
snri = (-10:1:28);
snrLin = 10.^(snri/10);
pam = [2,4,8]; % [2,4,8,16,32,64];
dX = 0.001; % for numeric calculations

%% Initiate variables

capacity_PAM_FiniteDPC  = zeros(length(pam),length(snri));
capacity_UNI            = zeros(1,length(snri));
capacity_UNI_DPC        = zeros(1,length(snri));
capacity_UNI_FiniteDPC  = zeros(1,length(snri));
capacity_UNI_FiniteDPC_upperBound = zeros(1,length(snri));

%% Start loop for PAM
k = 3;
% Generate PAM augmented constellation outside the main for loop
for i=1:length(pam)
    Delta = pam(i) + dX;
    xAxis = (-(2*Delta)*k - Delta):dX : ((2*Delta)*k + Delta);
    xAxisForUniformDist = -Delta:dX:Delta;
    xAxisForModulo = 2*(-(2*Delta)*k - Delta):dX : 2*((2*Delta)*k + Delta);
    
    % indeces for modulo
    if i==1
        idxArray = zeros(4*k+1,length(xAxisForUniformDist)-1);
        idxArrayNoise = zeros(2*k+1,length(xAxisForUniformDist)-1);
        cnt1 = 1; cnt2 = 1;
        for m = -1*2*k:2*k
            temp1 = find((xAxisForModulo >= (2*Delta*m - Delta)) & (xAxisForModulo <= (2*Delta*m + Delta)) & (xAxisForModulo ~= 0));
            idxArray(cnt1,:) = temp1(1,1:size(idxArray,2));
            if ((m >= -1*k) && (m <= 1*k))
                temp2 = find((xAxis >= (2*Delta*m - Delta)) & (xAxis <= (2*Delta*m + Delta)) & (xAxis ~= 0));
                idxArrayNoise(cnt2,:) = temp2(1,1:size(idxArrayNoise,2));
                cnt2 = cnt2 + 1;
            end
            cnt1 = cnt1 + 1;
        end
    end
    %% inputs - pam
    pamSymbols = (-1*(pam(i) - 1) : 2 : (pam(i) - 1));
    pamOrigDist = ones(size(pamSymbols))/pam(i);
    pamDist = zeros(size(xAxis));
    tempPamIdx = sum(abs(xAxis(:) - pamSymbols) < dX/2,2);
    pamDist(tempPamIdx > 0.5) = 1/pam(i);
    
    % augmented constellation : (v tilde)augmented pam with non uniform probabilities;
    augPamDist = zeros(size(xAxis));
    augPam = [((1 : 2 : (pam(i) - 1)) - 2*Delta) (-1*(pam(i) - 1) : 2 : (pam(i) - 1)) ((-1*(pam(i) - 1) : 2 : -1) + 2*Delta)];
    tempPamIdx = sum(abs(xAxis(:) - augPam) < dX/2,2);
    tempPamProbs = [ abs((-1*(pam(i) - 1) : 2 : -1))/(2*Delta) ...
        (2*Delta - abs((-1*(pam(i) - 1) : 2 : (pam(i) - 1))))/(2*Delta) ...
        abs((-1*(pam(i) - 1) : 2 : -1))/(2*Delta)]/pam(i);
    augPamDist(tempPamIdx > 0.5) = tempPamProbs/dX;
    
    % Calculate composite noise distribution
    pam_oneProbability = zeros(size(pamSymbols));pam_minusOneProbability = zeros(size(pamSymbols));
    
    pam_oneProbability(pamSymbols >= 0) = abs(pamSymbols(pamSymbols >= 0))/(2*Delta);
    pam_minusOneProbability(pamSymbols <= 0) = abs(pamSymbols(pamSymbols <= 0))/(2*Delta);
    
    pam_probs = [1-abs(pamSymbols(:))/(2*Delta) pam_oneProbability(:) pam_minusOneProbability(:)];
    pam_centers = [0*ones(length(pamSymbols),1) 2*Delta*ones(length(pamSymbols),1) -1*2*Delta*ones(length(pamSymbols),1)];
    
    % calculate original probabilities entropy
    pam_logProbs = log2(pam_probs); pam_logProbs(pam_logProbs == -inf) = 0;
    pam_origConditionedEntropy = sum(pam_logProbs.*pam_probs,2);
    pam_naiveNoiseEntropy = -1*dX*(pamOrigDist * pam_origConditionedEntropy);
    
    %% inputs - uniform
    % uniform
    if (i==1)
        uniformDist = zeros(size(xAxis));
        uniformDist(xAxis >= -1*Delta & xAxis <= Delta) = 1/(2*Delta);
        
        % augmented constellation
        triangleDist = conv(uniformDist,uniformDist)*dX;
        ValidIdx = round(length(uniformDist)/2):(length(triangleDist)-round(length(uniformDist)/2));
        triangleDist = triangleDist(ValidIdx);
        
        % Calculate composite noise distribution
        uni_oneProbability = zeros(size(xAxisForUniformDist));uni_minusOneProbability = zeros(size(xAxisForUniformDist));
        
        uni_oneProbability(xAxisForUniformDist >= 0) = abs(xAxisForUniformDist(xAxisForUniformDist >= 0))/(2*Delta);
        uni_minusOneProbability(xAxisForUniformDist <= 0) = abs(xAxisForUniformDist(xAxisForUniformDist <= 0))/(2*Delta);
        
        uni_probs = [1-abs(xAxisForUniformDist(:))/(2*Delta) uni_oneProbability(:) uni_minusOneProbability(:)];
        uni_centers = [0*ones(length(xAxisForUniformDist),1) 2*Delta*ones(length(xAxisForUniformDist),1)...
            -1*2*Delta*ones(length(xAxisForUniformDist),1)];
        
        % calculate entropy using gauss-hermite approx. set minimum to original
        % probabilities entropy
        uni_logProbs = log2(uni_probs); uni_logProbs(uni_logProbs == -inf) = 0;
        origConditionedEntropy = sum(uni_logProbs.*uni_probs,2);
        uni_naiveNoiseEntropy = -1*dX*sum(origConditionedEntropy)/(2*Delta);
        
    end
    
    % loop over SNRs
    for j=1:length(snri)
        sigma = sqrt((Delta^2/3)*10^(-snri(j)/10));
        
        %% calculate distributions
        % noise
        noiseDist = (1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/(sigma^2));
        
        % calculate channel output distribution
        if i==1
            % calculate equivalent noise for the DPC case alpha*X + (1-alpha)*N where alpha = (SNR/(SNR+1))
            alpha = snrLin(j)/(1+snrLin(j));
            alphaTimesNoise = (1/alpha)*(1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/((alpha^2)*(sigma^2)));
            alphaTimesInput = zeros(size(xAxis));
            alphaTimesInput(xAxis >= -1*Delta*(1-alpha) & xAxis <= Delta*(1-alpha)) = 1/(2*Delta*(1-alpha));
            
            equivNoiseDPC = conv(alphaTimesNoise,alphaTimesInput) * dX;
            
            
            % regulat uniform input in awgn channel
            yDist = conv(noiseDist,uniformDist) * dX;
            
            % uniform input + uniform interfere , modulo in Tx
            yDistFiniteDPC_uni = conv(noiseDist,triangleDist) * dX;
            
            % DPC : uniform input + noise and modulo
            yDistDPC = sum(yDist(idxArray),1);
            noiseDistDPC = sum(equivNoiseDPC(idxArray),1);
        end
        yDistFiniteDPC_pam = conv(noiseDist,augPamDist) * dX;
        
        %% calculate H(Y|X) - first condition on X=x and then average
        % calculation for the finite case
        uni_noiseConditionedEntropy = gaussianMixtureEntropyApprox(uni_centers,uni_probs,sigma,20);
        pam_noiseConditionedEntropy = gaussianMixtureEntropyApprox(pam_centers,pam_probs,sigma,20);
        
        uni_curNoiseEntropyFiniteDPC = -1*dX*sum(uni_noiseConditionedEntropy)/(2*Delta);
        pam_curNoiseEntropyFiniteDPC = -1*sum(pam_noiseConditionedEntropy)/(pam(i));
        
        % calculate for the regular DPC case
        if i==1
            curNoiseEntropyDPC = -1*sum(noiseDistDPC(noiseDistDPC > 0).*log2(noiseDistDPC(noiseDistDPC > 0))) * dX;
        end
        
        %% calculate H(Y) : numeric integration
        curEntropyAWGN_uni      = -1*sum(yDist(yDist > 0).*log2(yDist(yDist > 0))) * dX;
        curEntropyDPC_uni       = -1*sum(yDistDPC(yDistDPC > 0).*log2(yDistDPC(yDistDPC > 0))) * dX;
        curEntropyFiniteDPC_uni = -1*sum(yDistFiniteDPC_uni(yDistFiniteDPC_uni > 0).*log2(yDistFiniteDPC_uni(yDistFiniteDPC_uni > 0))) * dX;
        curEntropyFiniteDPC_pam = -1*sum(yDistFiniteDPC_pam(yDistFiniteDPC_pam > 0).*log2(yDistFiniteDPC_pam(yDistFiniteDPC_pam > 0))) * dX;
        
        % calculate MI for current SNR
        if i==1
            capacity_UNI(j)   = (curEntropyAWGN_uni - 0.5*log2(2*pi*exp(1)*sigma^2));
            capacity_UNI_DPC(j)   = (curEntropyDPC_uni - curNoiseEntropyDPC);
            
            if isnan(uni_curNoiseEntropyFiniteDPC)
                capacity_UNI_FiniteDPC(j) = (curEntropyFiniteDPC_uni - (uni_naiveNoiseEntropy + 0.5*log2(2*pi*exp(1)*sigma^2)));
            else
                capacity_UNI_FiniteDPC(j) = (curEntropyFiniteDPC_uni - uni_curNoiseEntropyFiniteDPC);
            end
            dt = 1e-4;
            t = dt:dt:0.5;
            capacity_UNI_FiniteDPC_upperBound(j) = 0.5*log2(1+2*snrLin(j)) + 0.223 + ...
            2*dt*sum(t.*log2(t*(1-exp(-1/(4*sigma^2))) + exp(-1/(4*sigma^2))) + ...
            (1-t).*log2(1 - t*(1-exp(-1/(4*sigma^2)))));
        end
        
        if isnan(pam_curNoiseEntropyFiniteDPC)
            capacity_PAM_FiniteDPC(i,j) = log2(pam(i));
        else
            capacity_PAM_FiniteDPC(i,j) = (curEntropyFiniteDPC_pam - pam_curNoiseEntropyFiniteDPC);
        end
        
        
    end
    display(strcat('Finished PAM : ',num2str(pam(i))));
end

% clip at the minimum to avoid negative / zeros (that arise from numeric shit)
capacity_UNI = max(capacity_UNI,min(abs(capacity_UNI)));
capacity_UNI_DPC = max(capacity_UNI_DPC,min(abs(capacity_UNI_DPC)));
capacity_UNI_FiniteDPC = max(capacity_UNI_FiniteDPC,min(min(abs(capacity_UNI_FiniteDPC))));
capacity_PAM_FiniteDPC = capacity_PAM_FiniteDPC.*sign(capacity_PAM_FiniteDPC);
capacity_UNI_FiniteDPC_upperBound = capacity_UNI_FiniteDPC_upperBound.*sign(capacity_UNI_FiniteDPC_upperBound);

% AWGN Channel capacity
GaussianCapacity = 0.5*log2(1+ 10.^(snri/10));
lowerBound = max(0,0.5*log2(1+ 2 * 10.^(snri/10)) - 0.724);
upperBound = capacity_UNI_FiniteDPC_upperBound;

%% Plot Graphs - as function of SNR and Eb/N0
% Bounds and awgn Capacity
figure;hold all
plot(snri,GaussianCapacity,'-*','LineWidth',1.6);
plot(snri,capacity_UNI,'-o','LineWidth',1.6);
plot(snri,capacity_UNI_FiniteDPC,'-p','LineWidth',1.6);
plot(snri,capacity_UNI_DPC,'-^','LineWidth',1.6);
plot(snri,lowerBound,'--','LineWidth',1.6);
plot(snri,upperBound,'--','LineWidth',1.6);
set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('SNR (Es/N0) dB','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');

legend('Gaussian Capacity', 'Gaussian Capacity Uniform Input', ...
    'Finite Power DPC uniform interferer and input I(v;Y)','Regular DPC - Uniform Input',...
    'Finite DPC - Lower Bound','Finite DPC - Upper Bound');
title('Capacity for AWGN , With uniform interferer');

% plot Vs EsN0 - Uniform input and PAM input
figure;hold all
plot(snri,GaussianCapacity,'-*','LineWidth',1.6);
plot(snri,capacity_UNI,'-o','LineWidth',1.6);
plot(snri,capacity_UNI_FiniteDPC,'-p','LineWidth',1.6);
plot(snri,capacity_UNI_DPC,'-^','LineWidth',1.6);
for i=1:length(pam)
    plot(snri,capacity_PAM_FiniteDPC(i,:),'LineWidth',1.6);
end
set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('SNR (Es/N0) dB','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');

legend('Gaussian Capacity', 'Gaussian Capacity Uniform Input', ...
    'Finite Power DPC uniform interferer and input I(v;Y)','Regular DPC - Uniform Input',...
    'Finite DPC : 2PAM','Finite DPC : 4PAM','Finite DPC : 8PAM');
title('Capacity for AWGN , With uniform interferer');

% plot Vs EbN0
figure;hold all
plot(10*log10(snrLin./(2*GaussianCapacity)),GaussianCapacity,'-*','LineWidth',1.6);
plot(10*log10(snrLin./(2*capacity_UNI)),capacity_UNI,'-o','LineWidth',1.6);
plot(10*log10(snrLin./(2*capacity_UNI_FiniteDPC)),capacity_UNI_FiniteDPC,'--','LineWidth',1.6);
plot(10*log10(snrLin./(2*capacity_UNI_DPC)),capacity_UNI_DPC,'LineWidth',1.6);
for i=1:length(pam)
    plot(10*log10(snrLin./(2*capacity_PAM_FiniteDPC(i,:))),capacity_PAM_FiniteDPC(i,:),'LineWidth',1.6);
end
set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('SNR (Eb/N0) dB','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');

legend('Gaussian Capacity', 'Gaussian Capacity Uniform Input', ...
    'Finite Power DPC uniform interferer and input I(v;Y)','Regular DPC - Uniform Input',...
    '2PAM','4PAM','8PAM');
title('Capacity for AWGN , With uniform interferer');
