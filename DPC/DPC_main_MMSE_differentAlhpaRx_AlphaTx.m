addpath(genpath('..\Utils'))

% This script calculate capacity for different DPC settings. The main
% purpose is to calculate the capacity with finite power interferer ,
% without modulo in the reciever.
clc; clear ;
%% data
snri = -7.5;
snrLin = 10.^(snri/10);
inputType = {'uni'};
dX = 5e-4; % for numeric calculations

%% Initiate variables

% finite dpc (no modulo in the reciever)
capacity_FiniteDPC  = zeros(length(inputType),length(snri));
% regula dpc
capacity_DPC        = zeros(length(inputType),length(snri));
% awgn
capacity_AWGN       = zeros(length(inputType),length(snri));

alphaRx = 0.05:0.1:1;
alphaTx = 0.05:0.1:1;
%% Iterate over input
numOfSegments = 3;
for i=1:length(inputType)
    
    % generate axis from input
    [Delta,xAxis,ModuloIndices,normFactor] = parseInput(inputType(i),numOfSegments,dX);
    
    % loop over SNRs
    % regular DPC :
    % Y = modulo(v + (alpha-1)*X + alpha*Z) where
    %   o v is the information signal
    %   o X is the Tx output , which is uniform (becuase of the dither)
    %   o Z is the noise
    % Capacity is H(Y) - H(modulo(Z_eq = (alpha-1)*X + alpha*Z)
    % optimal alpha is the MMSE estimation coefficient
    % alpha = snr/(1+snr)
    
    % finite DPC :
    %  Y = alpha_Rx * mod(v - alpha_Tx*S+U) - U + alpha_Rx*S + alpha_Rx*Z
    %   o v is the information signal
    %   o k is the integer part from the modulo subtraction and addition of S
    %   o Z is the noise
    %  the dither makes lot of things uniform and independant
    % Capacity is H(Y) - H(Y|v)
    % The calculate is direct from the definition of entropies and signals
    % statistics
    
    for j=1:length(alphaRx)
        for k=1:length(alphaTx)
            
            % calculate optimal alphas
            alpha = snrLin/(1+snrLin);
            alpha_Tx = alphaTx(k);
            alpha_Rx = alphaRx(j);
            
            
            uniformDist = zeros(size(xAxis));
            uniformDist(abs(xAxis) < Delta) = 1/(2*Delta);
            
            uniformAlphaDist = zeros(size(xAxis));
            uniformAlphaDist(abs(xAxis) < Delta*alpha) = 1/(2*Delta*alpha);
            
            uniformAlphaRxDist = zeros(size(xAxis));
            uniformAlphaRxDist(abs(xAxis) < Delta*alpha_Rx) = 1/(2*Delta*alpha_Rx);
            
            % snr is relative to Tx Output : For the uniform interferer case Tx Output is
            % unifrom between [-Delta,Delta] , so Tx Power = Delta^2 / 3
            sigma = sqrt((Delta^2/3)*10^(-snri/10));
            
            %% calculate distributions
            % AWGN with Zero Mean
            noiseDist = (1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/(sigma^2));
            alphaTimesNoise = (1/alpha)*(1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/((alpha^2)*(sigma^2)));
            alphaRxTimesNoise = (1/alpha_Rx)*(1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/((alpha_Rx^2)*(sigma^2)));
            
            % channel output distribution :
            
            % uniform input in awgn
            yDist = conv(noiseDist,uniformDist) * dX;
            yDist = yDist (floor(length(noiseDist)/2) : end - floor(length(noiseDist)/2));
            
            % Finite Power DPC : Modulo only in Tx
            calc1 = conv(uniformAlphaRxDist,alphaRxTimesNoise)*dX;
            calc1 = calc1(floor(length(uniformAlphaRxDist)/2) : end - floor(length(uniformAlphaRxDist)/2));
            temp_yDist_FiniteDPC = conv(calc1,uniformAlphaRxDist) * dX;
            temp_yDist_FiniteDPC = temp_yDist_FiniteDPC(floor(length(calc1)/2) : end - floor(length(calc1)/2));
            yDist_FiniteDPC = conv(temp_yDist_FiniteDPC,uniformDist) * dX;
            yDist_FiniteDPC = yDist_FiniteDPC(floor(length(temp_yDist_FiniteDPC)/2) : end - floor(length(temp_yDist_FiniteDPC)/2));
            
            alphaTimesInput = zeros(size(xAxis));
            alphaTimesInput(xAxis >= -1*Delta*(1-alpha) & xAxis <= Delta*(1-alpha)) = 1/(2*Delta*(1-alpha));
            equivNoiseDPC = conv(alphaTimesNoise,alphaTimesInput) * dX;
            equivNoiseDPC = equivNoiseDPC(ceil(length(alphaTimesNoise)/2) : end - ceil(length(alphaTimesNoise)/2));
            noiseDistDPC = sum(equivNoiseDPC(ModuloIndices),2);
            
            
            %% calculate H(Y|X) - first condition on X=x and then average
            % calculation for the finite case
            noiseConditionedEntropy = calcFiniteDPC_entropy_alphaRxTx(Delta,alpha_Tx,alpha_Rx,sigma);
            noiseEntropy_FiniteDPC = noiseConditionedEntropy;
            
            noiseEntropy_DPC = -1*sum(noiseDistDPC(noiseDistDPC > 0).*log2(noiseDistDPC(noiseDistDPC > 0))) * dX;
            
            
            %% calculate H(Y) : numeric integration
            curEntropy_DPC          = -1*sum(uniformDist(uniformDist > 0).*log2(uniformDist(uniformDist > 0))) * dX;
            curEntropy_FiniteDPC    = -1*sum(yDist_FiniteDPC(yDist_FiniteDPC > 0).*log2(yDist_FiniteDPC(yDist_FiniteDPC > 0))) * dX;
            
            % calculate MI for current SNR
            capacity_FiniteDPC(i,j,k) = curEntropy_FiniteDPC - noiseEntropy_FiniteDPC;
            capacity_DPC(i,j,k)       = curEntropy_DPC - noiseEntropy_DPC;
            % calculate the regular uniform input in awgn capacity
            curEntropyAWGN_uni      = -1*sum(yDist(yDist > 0).*log2(yDist(yDist > 0))) * dX;
            capacity_AWGN(i,j,k) = curEntropyAWGN_uni - 0.5*log2(2*pi*exp(1)*sigma^2);
        end
    end
    display(strcat('Finished : ',inputType(i)));
    
end

% clip at the minimum to avoid negative / zeros (that arise from numeric shit)
[capacity_AWGN,capacity_DPC,capacity_FiniteDPC] = clipCapacity(capacity_AWGN,capacity_DPC,capacity_FiniteDPC);

resultPlot(snri,capacity_AWGN,capacity_DPC,capacity_FiniteDPC,inputType,alphaRx,alphaTx);


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

function [capacity_AWGN,capacity_DPC,capacity_FiniteDPC] = clipCapacity(capacity_AWGN,capacity_DPC,capacity_FiniteDPC)

capacity_AWGN = reshape(max(capacity_AWGN(:),abs(min(capacity_AWGN(:)))),size(capacity_AWGN));
capacity_DPC = reshape(max(capacity_DPC(:),abs(min(capacity_DPC(:)))),size(capacity_DPC));
capacity_FiniteDPC = reshape(max(capacity_FiniteDPC(:),abs(min(capacity_FiniteDPC(:)))),size(capacity_FiniteDPC));


end

function [] =   resultPlot(snr,capacity_AWGN,capacity_DPC,capacity_FiniteDPC,inputType,alphaRx,alphaTx)

snrLin = 10.^(snr/10);
% AWGN Channel capacity
GaussianCapacity = 0.5*log2(1+ 10.^(snr/10));

[alphaRx_Plot,alphaTx_Plot] = meshgrid(alphaRx,alphaTx);
% Plot Finite DPC + DPC + AWGN + Lower Bound , for uniform input , EsN0
uniform_idx = find(strcmp(inputType,'uni'));
% Bounds and awgn Capacity
figure;hold all
plot3(alphaRx_Plot,alphaTx_Plot,GaussianCapacity*ones(size(alphaRx)),'-*','LineWidth',1.6);
plot3(alphaRx_Plot,alphaTx_Plot,capacity_AWGN(uniform_idx,:,:),'-o','LineWidth',1.6);
plot3(alphaRx_Plot,alphaTx_Plot,capacity_DPC(uniform_idx,:,:),'-p','LineWidth',1.6);
plot3(alphaRx_Plot,alphaTx_Plot,capacity_FiniteDPC(uniform_idx,:,:),'-^','LineWidth',1.6);

set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('\alpha_{Rx}','FontSize', 14, 'FontName', 'Times New Roman');
xlabel('\alpha_{Tx}','FontSize', 14, 'FontName', 'Times New Roman');
zlabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');

legend('Gaussian Capacity', 'Gaussian Capacity Uniform Input', ...
    'DPC uniform interferer and input I(v;Y)','Finite Power DPC - Uniform Input');
title('Capacity for AWGN , With uniform interferer and different \alpha_{Rx} , \alpha_{Tx}');

end

