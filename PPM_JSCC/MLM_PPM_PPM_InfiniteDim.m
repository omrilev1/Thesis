%% Simulation of Successive MLM and Linear transmissions, for the transmissions of an uniform source
%% over a channel with unlimited bandwidth. We optimize the scheme parameters according to the simulated profile

close all; clear all; clc;

profileOrder = 2;
saveResults = 1;

% Init parameters and arrays structures
SNR = -7:0.25:20;  % 8              % Direct channel SNR
snrLin = 10.^(SNR/10);
profile = 1./(1 + (snrLin).^profileOrder);

% profile = 1./(1 + (snrLin/snrLin(1)).^profileOrder);

maxStages = 20; % number of PPM levels

%% Modulo parameters:
% interval
% Delta = 0.88; alpha = 0.51;
% M = 1;
% energyAlloc = ones(1,M + 1) .* Delta .* exp(-(0:1:(M)) * alpha);
% energyAlloc(end) = 1.35*energyAlloc(end);

%
M = 1;
energyAlloc = [0.85 0.7];
disp(strcat('Total Energy = ',num2str(sum(energyAlloc))))

%% Simulation parameters

D_PPM = zeros(length(SNR),1);
SDR_tot = ones(length(SNR),1);
totalEnergy = ones(length(SNR),1);
currNumOfLevels = 1;
% currE = zeros(1,currNumOfLevels);currE(1) = 1;

etaLinear = 1;

for i=1:length(SNR)
    
    prevNumOfLevels = currNumOfLevels;
    if ((i > 1) && (currNumOfLevels <= M))
        if ((D_PPM(i - 1) > profile(i - 1)))
            currNumOfLevels = prevNumOfLevels + 1;
        end
    end
    
    if prevNumOfLevels ~= currNumOfLevels
        startDist = D_PPM(i - 1);
        startENR = 2*energyAlloc(currNumOfLevels)*snrLin(i);
        if currNumOfLevels > M
            beta = 0.55;
            disp(strcat('beta opt = ',num2str(beta)));
            [startSDR,CorrCoeff,CurrPower] = simulateGaussianPPM(beta,2*energyAlloc(currNumOfLevels)*snrLin(i));
            alpha_c0 = (CorrCoeff)/CurrPower;
        end
    end
    
    if currNumOfLevels == 1
        D_PPM(i) = 1/(1 + 2*energyAlloc(1)*snrLin(i));
    else
        
        if currNumOfLevels > M
            % simulate the SDR of PPM for the current SNR
            [SDR_PPM,CorrCoeff,CurrPower] = simulateGaussianPPM(beta,2*energyAlloc(currNumOfLevels)*snrLin(i));
            if SDR_PPM < 1
                SDR_PPM = 1;
                D_PPM(i) = (startDist/SDR_PPM);
            else
                noiseVar = 1 - (CorrCoeff)^2/CurrPower;
                
                %                 D_PPM(i) = (startDist*noiseVar)/(alpha_c0 +  (startDist/D_PPM(i-1)) * noiseVar);
                
                %                 D_PPM(i) = (startDist/(1 + SDR_PPM)) ...
                %                     * 1/(startSDR/(1 + startSDR) + (startDist/D_PPM(i-1)) * 1/(1 + SDR_PPM));
                
                %             D_PPM(i) = (startDist/SDR_PPM) ...
                %                 * 1/(1 + (startDist/D_PPM(i-1)) * 1/SDR_PPM);
                
                D_PPM(i) = (startDist/SDR_PPM);
            end
            SDR_tot(i) = SDR_PPM;
        else
            D_PPM(i) = ((1 + startENR)/startENR) * (startDist/(1 + 2*energyAlloc(currNumOfLevels)*snrLin(i)));
        end
    end
    totalEnergy(i) = sum(energyAlloc(1:currNumOfLevels));
    
    if mod(i,5) == 0
        disp(strcat('Finished SNR = ',num2str(SNR(i))));
    end
end

if saveResults
    save(strcat('PpmInfDim_Profile',num2str(profileOrder)),'D_PPM','SNR','totalEnergy');
end

figure;
semilogy(SNR,D_PPM,'LineWidth',2.5); hold on;
semilogy(SNR,profile,'-.','LineWidth',2.5);
xlabel('ENR [dB]','FontSize',14); ylabel('Distortion [dB]','FontSize',14);
lgd = legend({'PPM + MLM','Profile'},'FontSize',16,'TextColor','Black');
grid on; grid minor;
title(strcat('Distortion and Profile, order = ',num2str(profileOrder)),'FontSize',14);

function [SDR_PPM,CorrCoeff,PowerCoeff] = simulateGaussianPPM(beta,SNRlin)

dt = 1/(325*beta);

overload = 6.75;
t = -overload:dt:overload;

Nrun =  1e4;

ppmPulse = zeros(size(t));
ppmPulse(abs(t) < 1/(2*beta)) = sqrt(beta);
ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));
% ppmPulse(abs(t) > 0.5) = 0;

Lx = length(ppmPulse);
Ly = length(ppmPulse)+length(ppmPulse)-1;
Ly2 = pow2(nextpow2(Ly));
PPMfreq = fft(ppmPulse, Ly2);
currMSE_MAP = zeros(1,Nrun);
currCorr_MAP = zeros(1,Nrun);
currPower_MAP = zeros(1,Nrun);

parfor n=1:Nrun
    
    % generate source - Gaussian Source with edge truncation
    S = randn;
    if abs(S) > max(t)
        S = max(t)*sign(S);
    end
    
    TxPulse = zeros(size(t));
    TxPulse(abs(t - S) < 1/(2*beta)) = sqrt(beta);
    TxPulse = sqrt(SNRlin)*TxPulse/sum(abs(TxPulse.^2)*dt);
    
    noise = randn(size(t));
    noise = sqrt(1/(2*dt))*noise;
    r = TxPulse + noise;
    
    % PPM Correlator receiver
    PPMcorr = fconv(r,Lx,Ly,Ly2,PPMfreq);
    
    [~,maxIdx] = max(sqrt(SNRlin)*PPMcorr*dt - 0.5*(t.^2));
    sHat_MAP = t(maxIdx);
    
    currMSE_MAP(n) = (S - sHat_MAP)^2;
    currCorr_MAP(n) = S*(sHat_MAP);
    currPower_MAP(n) = (sHat_MAP)^2;
end

SDR_PPM = 1 / (sum(currMSE_MAP)/Nrun);
CorrCoeff = sum(currCorr_MAP)/Nrun;
PowerCoeff = sum(currPower_MAP)/Nrun;

end

function y = fconv(x,Lx,Ly,Ly2,H)
% Convolution in frequency domain using power of 2 fft
% Since the input signal is real we use known fft identities to accelerate
% the fft

% input fft
X = fft(x, Ly2);

% multiply with precalculated freq domain signal
Y = X.*H;

% inverse fft and truncation
y = real(ifft(Y, Ly2));
y = y(1:1:Ly);

if mod(Lx,2) == 0
    y = y(Lx/2 + 1 : end - (Lx/2) + 1);
else
    y = y((Lx+1)/2 : end - ((Lx+1)/2) + 1);
end

end

