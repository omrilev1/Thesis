addpath(genpath('..\Utils'))

% This script calculate capacity for different DPC settings. The main
% purpose is to calculate the capacity with finite power interferer ,
% without modulo in the reciever.
clc; clear ; close all;
%% data
snri = (-10:1:20);
snrLin = 10.^(snri/10);
SIR = 6; % SIR _ dB
SIR_Lin = 10^(-SIR/20);
inputType = {'uni'};
dX = 5e-4; % for numeric calculations

%% Initiate variables

% finite dpc (no modulo in the reciever)
capacity_FiniteDPC  = zeros(length(inputType),length(snri));
% regula dpc
capacity_DPC        = zeros(length(inputType),length(snri));
% awgn
capacity_AWGN       = zeros(length(inputType),length(snri));


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
    
    % finite DPC :
    %  Y = alpha * mod(v - alpha*S+U) - U + alpha*S + alpha*Z
    %   o v is the information signal
    %   o k is the integer part from the modulo subtraction and addition of S
    %   o Z is the noise
    %  the dither makes lot of things uniform and independant
    % Capacity is H(Y) - H(Y|v)
    % The calculate is direct from the definition of entropies and signals
    % statistics
    
    for j=1:length(snri)
        
        alpha = snrLin(j)/(1+snrLin(j));
        
        uniformDist = zeros(size(xAxis));
        uniformDist(abs(xAxis) < Delta) = 1/(2*Delta);
        
        uniformAlphaDist = zeros(size(xAxis));
        uniformAlphaDist(abs(xAxis) < Delta*alpha) = 1/(2*Delta*alpha);
        
        % generate interferer dist 
        uniformAlpha_IntrfrDist = zeros(size(xAxis));
        uniformAlpha_IntrfrDist(abs(xAxis) < Delta*alpha*SIR_Lin) = 1/(2*Delta*alpha*SIR_Lin);        
        % snr is relative to Tx Output : For the uniform interferer case Tx Output is
        % unifrom between [-Delta,Delta] , so Tx Power = Delta^2 / 3
        sigma = sqrt((Delta^2/3)*10^(-snri(j)/10));
        
        %% calculate distributions
        % AWGN with Zero Mean
        noiseDist = (1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/(sigma^2));
        alphaTimesNoise = (1/alpha)*(1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/((alpha^2)*(sigma^2)));
        
        % channel output distribution :
        
        % uniform input in awgn
        yDist = conv(noiseDist,uniformDist) * dX;
        yDist = yDist (floor(length(noiseDist)/2) : end - floor(length(noiseDist)/2));
        
        % Finite Power DPC : Modulo only in Tx
        calc1 = conv(uniformAlphaDist,alphaTimesNoise)*dX;
        calc1 = calc1(floor(length(uniformAlphaDist)/2) : end - floor(length(uniformAlphaDist)/2));
        temp_yDist_FiniteDPC = conv(calc1,uniformAlpha_IntrfrDist) * dX;
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
        noiseConditionedEntropy = calcFiniteDPC_entropy(Delta,alpha,sigma,SIR);
        noiseEntropy_FiniteDPC = noiseConditionedEntropy;
        
        noiseEntropy_DPC = -1*sum(noiseDistDPC(noiseDistDPC > 0).*log2(noiseDistDPC(noiseDistDPC > 0))) * dX;
        
        
        %% calculate H(Y) : numeric integration
        curEntropy_DPC          = -1*sum(uniformDist(uniformDist > 0).*log2(uniformDist(uniformDist > 0))) * dX;
        curEntropy_FiniteDPC    = -1*sum(yDist_FiniteDPC(yDist_FiniteDPC > 0).*log2(yDist_FiniteDPC(yDist_FiniteDPC > 0))) * dX;
        
        % calculate MI for current SNR
        capacity_FiniteDPC(i,j) = curEntropy_FiniteDPC - noiseEntropy_FiniteDPC;
        capacity_DPC(i,j)       = curEntropy_DPC - noiseEntropy_DPC;
        % calculate the regular uniform input in awgn capacity
        curEntropyAWGN_uni      = -1*sum(yDist(yDist > 0).*log2(yDist(yDist > 0))) * dX;
        capacity_AWGN(i,j) = curEntropyAWGN_uni - 0.5*log2(2*pi*exp(1)*sigma^2);
    end
    
    display(strcat('Finished : ',inputType(i)));
    
end

% clip at the minimum to avoid negative / zeros (that arise from numeric shit)
[capacity_AWGN,capacity_DPC,capacity_FiniteDPC] = clipCapacity(capacity_AWGN,capacity_DPC,capacity_FiniteDPC);

resultPlot(snri,capacity_AWGN,capacity_DPC,capacity_FiniteDPC,inputType,SIR);


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

function [] =   resultPlot(snr,capacity_AWGN,capacity_DPC,capacity_FiniteDPC,inputType,SIR)

snrLin = 10.^(snr/10);
% AWGN Channel capacity
GaussianCapacity = 0.5*log2(1+ 10.^(snr/10));

% Plot Finite DPC + DPC + AWGN + Lower Bound , for uniform input , EsN0
uniform_idx = find(strcmp(inputType,'uni'));
% Bounds and awgn Capacity
figure;hold all
plot(snr,GaussianCapacity,'-*','LineWidth',2.5);
plot(snr,capacity_AWGN(uniform_idx,:),'-o','LineWidth',2.5);
plot(snr,capacity_DPC(uniform_idx,:),'-p','LineWidth',2.5);
plot(snr,capacity_FiniteDPC(uniform_idx,:),'-^','LineWidth',2.5);

set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('SNR (Es/N0) [dB]','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');

legend('Gaussian Capacity', 'Gaussian Capacity Uniform Input', ...
    'DPC uniform interferer and input I(v;Y)','Finite Power DPC - Uniform Input');
title({'Capacity for AWGN , With uniform interferer and alpha_{MMSE}',strcat('SIR = ',num2str(SIR),' [dB]')});

% Plot Vs EbN0
% Bounds and awgn Capacity
figure;hold all
plot(10*log10(snrLin./(2*GaussianCapacity)),GaussianCapacity,'-*','LineWidth',2.5);
plot(10*log10(snrLin./(2*capacity_AWGN(uniform_idx,:))),capacity_AWGN(uniform_idx,:),'-o','LineWidth',2.5);
plot(10*log10(snrLin./(2*capacity_DPC(uniform_idx,:))),capacity_DPC(uniform_idx,:),'-p','LineWidth',2.5);
plot(10*log10(snrLin./(2*capacity_FiniteDPC(uniform_idx,:))),capacity_FiniteDPC(uniform_idx,:),'-^','LineWidth',2.5);

set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('Eb / N0 [dB]','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');

legend('Gaussian Capacity', 'Gaussian Capacity Uniform Input', ...
    'DPC uniform interferer and input I(v;Y)','Finite Power DPC - Uniform Input');
title({'Capacity for AWGN , With uniform interferer and alpha_{MMSE}',strcat('SIR = ',num2str(SIR),' [dB]')});
end

