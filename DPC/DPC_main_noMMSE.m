addpath(genpath('..\Utils'))

% This script calculate capacity for different DPC settings. The main
% purpose is to calculate the capacity with finite power interferer ,
% without modulo in the reciever.
clc; clear ; close all;
%% data
snri = (-10:1:28);
snrLin = 10.^(snri/10);
inputType = {'uni'}; % [2,4,8,16,32,64];
toUseMMSE = [0];

dX = 1e-4; % for numeric calculations

%% Initiate variables

% finite dpc (no modulo in the reciever)
capacity_FiniteDPC  = zeros(length(inputType),length(snri),length(toUseMMSE));
% regula dpc
capacity_DPC        = zeros(length(inputType),length(snri),length(toUseMMSE));
% awgn
capacity_AWGN       = zeros(length(inputType),length(snri));


%% Iterate over input
numOfSegments = 3;
for i=1:length(inputType)
    
    % generate axis from input
    [Delta,xAxis,ModuloIndices,normFactor] = parseInput(inputType(i),numOfSegments,dX);
    
    % calculate input distribution , and augmented constellation distribution (v tilde)
    [origPDF,augPDF,intPart_Val,intPart_prob,naiveNoiseEntropy] = calcConstellations(inputType(i),Delta,xAxis,dX);
    
    uniformDist = zeros(size(xAxis));
    uniformDist(abs(xAxis) < Delta) = 1/(2*Delta);
    % loop over SNRs
    for j=1:length(snri)
        
        % snr is relative to Tx Output : For the uniform interferer case Tx Output is
        % unifrom between [-Delta,Delta] , so Tx Power = Delta^2 / 3
        sigma = sqrt((Delta^2/3)*10^(-snri(j)/10));
        
        %% calculate distributions
        % AWGN with Zero Mean
        noiseDist = (1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/(sigma^2));
        
        % channel output distribution :
        
        % uniform input in awgn
        yDist = conv(noiseDist,uniformDist) * dX;
        
        % Finite Power DPC : Modulo only in Tx
        yDist_FiniteDPC = conv(noiseDist,augPDF) * dX;
        
        % calculate with optimal MMSE Estimation Coefficient and without
        parfor k = 1 : length(toUseMMSE)
            % regular DPC :
            % Y = modulo(v + (alpha-1)*X + alpha*Z) where
            %   o v is the information signal
            %   o X is the Tx output , which is uniform (becuase of the dither)
            %   o Z is the noise
            % Capacity is H(Y) - H(modulo(Z_eq = (alpha-1)*X + alpha*Z)
            
            % finite DPC :
            %  Y = (modulo(v-S) + S + Z) = v_tilde + Z = v + 2Delta*k + Z where
            %   o v is the information signal
            %   o k is the integer part from the modulo subtraction and addition of S
            %   o Z is the noise
            % Capacity is H(Y) - H(2Delta*alpha*k + Z|v) = H(v_tilde + Z) - H(2Delta*alpha*k + Z|v)
            % in case we apply BLE in the reciever , we'll need to add log(alpha)/alpha for the capacity
            
            if toUseMMSE(k)
                alphaDPC = snrLin(j)/(1+snrLin(j));
                alpha_FiniteDPC = 2*snrLin(j)/(1+2*snrLin(j));
                
                alphaTimesNoise = (1/alphaDPC)*(1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/((alphaDPC^2)*(sigma^2)));
                alphaTimesInput = zeros(size(xAxis));
                alphaTimesInput(xAxis >= -1*Delta*(1-alphaDPC) & xAxis <= Delta*(1-alphaDPC)) = 1/(2*Delta*(1-alphaDPC));
                equivNoiseDPC = conv(alphaTimesNoise,alphaTimesInput) * dX;
                equivNoiseDPC = equivNoiseDPC(ceil(length(alphaTimesNoise)/2) : end - ceil(length(alphaTimesNoise)/2));
                noiseDistDPC = sum(equivNoiseDPC(ModuloIndices),2);
            else
                alpha_FiniteDPC = 1;
                alphaDPC = 1;
                noiseDistDPC = sum(noiseDist(ModuloIndices),2);
            end
            
            %% calculate H(Y|X) - first condition on X=x and then average
            % calculation for the finite case
            noiseConditionedEntropy = gaussianMixtureEntropyApprox(intPart_Val*alpha_FiniteDPC,intPart_prob,sigma*alpha_FiniteDPC,20);
            noiseEntropy_FiniteDPC = -1*sum(noiseConditionedEntropy)/normFactor;
            
            noiseEntropy_DPC = -1*sum(noiseDistDPC(noiseDistDPC > 0).*log2(noiseDistDPC(noiseDistDPC > 0))) * dX;
            
            
            %% calculate H(Y) : numeric integration
            curEntropy_DPC          = -1*sum(uniformDist(uniformDist > 0).*log2(uniformDist(uniformDist > 0))) * dX;
            curEntropy_FiniteDPC    = log2(abs(alpha_FiniteDPC)) ...
                -1*sum(yDist_FiniteDPC(yDist_FiniteDPC > 0).*log2(yDist_FiniteDPC(yDist_FiniteDPC > 0))) * dX;
            
            % calculate MI for current SNR
            if isnan(noiseEntropy_FiniteDPC)
                capacity_FiniteDPC(i,j,k) = curEntropy_FiniteDPC - naiveNoiseEntropy - 0.5*log2(2*pi*exp(1)*(alpha_FiniteDPC*sigma)^2);
            else
                capacity_FiniteDPC(i,j,k) = curEntropy_FiniteDPC - noiseEntropy_FiniteDPC;
            end
            
            capacity_DPC(i,j,k)       = curEntropy_DPC - noiseEntropy_DPC;
        end
        % calculate the regular uniform input in awgn capacity
        curEntropyAWGN_uni      = -1*sum(yDist(yDist > 0).*log2(yDist(yDist > 0))) * dX;
        capacity_AWGN(i,j) = curEntropyAWGN_uni - 0.5*log2(2*pi*exp(1)*sigma^2);
    end
    
    display(strcat('Finished : ',inputType(i)));
end

% clip at the minimum to avoid negative / zeros (that arise from numeric shit)
[capacity_AWGN,capacity_DPC,capacity_FiniteDPC] = clipCapacity(capacity_AWGN,capacity_DPC,capacity_FiniteDPC);

resultPlot(snri,capacity_AWGN,capacity_DPC,capacity_FiniteDPC,inputType,toUseMMSE);


%% Functions part
function [Delta,xAxis,ModuloIndices,normFactor] = parseInput(inputType,numOfSegments,dX)

charInput = char(inputType);
switch charInput(1:3)
    case 'uni'
        Delta = 1;
        normFactor = 2*Delta/dX;
    case 'PAM'
        Delta = eval(charInput(4:end)) - 1;
        normFactor = Delta + 1;
end
xAxis = (-(2*Delta)*numOfSegments - Delta): dX : ((2*Delta)*numOfSegments + Delta);
ModuloIndices = reshape(1:(length(xAxis) - 1),[],2*numOfSegments + 1);

end

function [origPDF,augPDF,intPart_Val,intPart_prob,intPart_Entropy] = calcConstellations(inputType,Delta,xAxis,dX)

charInput = char(inputType);
switch charInput(1:3)
    
    case 'uni'
        
        % original information signal distribution : unifrom over
        % [-Delta,Delta]
        origPDF = zeros(size(xAxis));
        origPDF(xAxis >= -1*Delta & xAxis <= Delta) = 1/(2*Delta);
        
        % augmented constellation : In scae of unifrom interferer we get
        % triangle between [-2*Delta,2*Delta]
        augPDF = conv(origPDF,origPDF)*dX;
        ValidIdx = round(length(origPDF)/2):(length(augPDF)-round(length(origPDF)/2) + 1);
        augPDF = augPDF(ValidIdx);
        
        % Calculate composite noise distribution
        AxisForUniform = -Delta:dX:Delta;
        probOne = zeros(size(AxisForUniform));probMinusOne = zeros(size(AxisForUniform));
        
        probOne(AxisForUniform >= 0) = abs(AxisForUniform(AxisForUniform >= 0))/(2*Delta);
        probMinusOne(AxisForUniform <= 0) = abs(AxisForUniform(AxisForUniform <= 0))/(2*Delta);
        
        intPart_prob = [1-abs(AxisForUniform(:))/(2*Delta) probOne(:) probMinusOne(:)];
        intPart_Val = [0*ones(length(AxisForUniform),1) 2*Delta*ones(length(AxisForUniform),1) -1*2*Delta*ones(length(AxisForUniform),1)];
        
        % calculate entropy using gauss-hermite approx. set minimum to original
        % probabilities entropy
        logProbs = log2(intPart_prob); logProbs(logProbs == -inf) = 0;
        conditionedEntropy = sum(logProbs.*intPart_prob,2);
        intPart_Entropy = -1*dX*sum(conditionedEntropy)/(2*Delta);
        
    case 'PAM'
        pamMax = Delta;
        %% inputs - pam : between [-Delta,Delta] with little margin to avoid numeric shit with modulo
        %  original probabilites are uniform over all PAM
        pamSymbols = -1*pamMax : 2 : pamMax ;
        pamSymbols(end) = pamSymbols(end) - dX;
        pamSymbols(1) = pamSymbols(1) + dX;
        
        pamOrigDist = ones(size(pamSymbols))/(Delta+1);
        origPDF = zeros(size(xAxis));
        tempPamIdx = sum(abs(xAxis(:) - pamSymbols) < dX/2,2);
        origPDF(tempPamIdx > 0.5) = 1/(Delta+1);
        
        % augmented constellation : (v tilde)augmented pam with non uniform probabilities;
        augPDF = zeros(size(xAxis));
        augPam = [((1 : 2 : pamMax) - 2*2*Delta) (-1*pamMax : 2 : pamMax) ((-1*pamMax : 2 : -1) + 2*2*Delta)];
        tempPamIdx = sum(abs(xAxis(:) - augPam) < dX/2,2);
        tempPamProbs = [ abs((-1*pamMax : 2 : -1))/(2*2*Delta) ...
            (2*2*Delta - abs((-1*pamMax : 2 : pamMax)))/(2*2*Delta) ...
            abs((-1*pamMax : 2 : -1))/(2*2*Delta)]/(2*Delta);
        augPDF(tempPamIdx > 0.5) = tempPamProbs/dX;
        
        % Calculate composite noise distribution
        probOne = zeros(size(pamSymbols));probMinusOne = zeros(size(pamSymbols));
        
        probOne(pamSymbols >= 0) = abs(pamSymbols(pamSymbols >= 0))/(2*2*Delta);
        probMinusOne(pamSymbols <= 0) = abs(pamSymbols(pamSymbols <= 0))/(2*2*Delta);
        
        intPart_prob = [1-abs(pamSymbols(:))/(2*2*Delta) probOne(:) probMinusOne(:)];
        intPart_Val = [0*ones(length(pamSymbols),1) 2*2*Delta*ones(length(pamSymbols),1) -1*2*2*Delta*ones(length(pamSymbols),1)];
        
        % calculate original probabilities entropy
        logProbs = log2(intPart_prob); logProbs(logProbs == -inf) = 0;
        conditionedEntropy = sum(logProbs.*intPart_prob,2);
        intPart_Entropy = -1*(pamOrigDist * conditionedEntropy);
        
end

end

function [capacity_AWGN,capacity_DPC,capacity_FiniteDPC] = clipCapacity(capacity_AWGN,capacity_DPC,capacity_FiniteDPC)

capacity_AWGN = reshape(max(capacity_AWGN(:),abs(min(capacity_AWGN(:)))),size(capacity_AWGN));
capacity_DPC = reshape(max(capacity_DPC(:),abs(min(capacity_DPC(:)))),size(capacity_DPC));
capacity_FiniteDPC = reshape(max(capacity_FiniteDPC(:),abs(min(capacity_FiniteDPC(:)))),size(capacity_FiniteDPC));


end

function [] =   resultPlot(snr,capacity_AWGN,capacity_DPC,capacity_FiniteDPC,inputType,toUseMMSE)

snrLin = 10.^(snr/10);
% AWGN Channel capacity
GaussianCapacity = 0.5*log2(1+ 10.^(snr/10));
lowerBound_FiniteDPC = max(0,0.5*log2(1+ 2 * 10.^(snr/10)) - 0.754);

% Plot Finite DPC + DPC + AWGN + Lower Bound , for uniform input , EsN0
uniform_idx = find(strcmp(inputType,'uni'));
for i = 1:length(toUseMMSE)
    % Bounds and awgn Capacity
    figure;hold all
    plot(snr,GaussianCapacity,'-*','LineWidth',1.6);
    plot(snr,capacity_AWGN(uniform_idx,:),'-o','LineWidth',1.6);
    plot(snr,capacity_DPC(uniform_idx,:,i),'-p','LineWidth',1.6);
    plot(snr,capacity_FiniteDPC(uniform_idx,:,i),'-^','LineWidth',1.6);
    plot(snr,lowerBound_FiniteDPC,'--','LineWidth',1.6);
    
    set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
    grid on; grid minor;
    box on;
    xlabel('SNR (Es/N0) [dB]','FontSize', 14, 'FontName', 'Times New Roman');
    ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');
    
    legend('Gaussian Capacity', 'Gaussian Capacity Uniform Input', ...
        'DPC uniform interferer and input I(v;Y)','Finite Power DPC - Uniform Input',...
        'Finite DPC - Lower Bound');
    title(strcat('Capacity for AWGN , With uniform interferer toUseMMSE = ',num2str(toUseMMSE(i))));
end
% Plot Vs EbN0
for i = 1:length(toUseMMSE)
    % Bounds and awgn Capacity
    figure;hold all
    plot(10*log10(snrLin./(2*GaussianCapacity)),GaussianCapacity,'-*','LineWidth',1.6);
    plot(10*log10(snrLin./(2*capacity_AWGN(uniform_idx,:))),capacity_AWGN(uniform_idx,:),'-o','LineWidth',1.6);
    plot(10*log10(snrLin./(2*capacity_DPC(uniform_idx,:,i))),capacity_DPC(uniform_idx,:,i),'-p','LineWidth',1.6);
    plot(10*log10(snrLin./(2*capacity_FiniteDPC(uniform_idx,:,i))),capacity_FiniteDPC(uniform_idx,:,i),'-^','LineWidth',1.6);
    plot(10*log10(snrLin./(2*lowerBound_FiniteDPC)),lowerBound_FiniteDPC,'--','LineWidth',1.6);
    
    set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
    grid on; grid minor;
    box on;
    xlabel('Eb / N0 [dB]','FontSize', 14, 'FontName', 'Times New Roman');
    ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');
    
    legend('Gaussian Capacity', 'Gaussian Capacity Uniform Input', ...
        'DPC uniform interferer and input I(v;Y)','Finite Power DPC - Uniform Input',...
        'Finite DPC - Lower Bound');
    title(strcat('Capacity for AWGN , With uniform interferer toUseMMSE = ',num2str(toUseMMSE(i))));
end

end

