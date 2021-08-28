%% Simulation of Successive MLM and Linear transmissions, for the transmissions of an uniform source
%% over a channel with unlimited bandwidth. We optimize the scheme parameters according to the simulated profile

close all; clear all; clc;

profileOrder = 2;
saveResults = 1;

% Init parameters and arrays structures
SNR = -7:0.25:30;  % 8              % Direct channel SNR
snrLin = 10.^(SNR/10);

profile = 1./(1 + (snrLin).^profileOrder);

maxStages = 50; % number of PPM levels

%% Modulo parameters:
% interval

Delta = 0.975; alpha = 0.65;
energyAlloc = ones(1,maxStages) .* Delta .* exp(-(0:1:(maxStages - 1)) * alpha);
disp(strcat('Total Energy = ',num2str(sum(energyAlloc))))

%% Simulation parameters

D_Linear = zeros(length(SNR),1);
totalEnergy = ones(length(SNR),1);

currNumOfLevels = 1;
% currE = zeros(1,currNumOfLevels);currE(1) = 1;

etaLinear = 1;

for i=1:length(SNR)
    
    prevNumOfLevels = currNumOfLevels;
    if i > 1
        if D_Linear(i - 1) >= profile(i - 1)
            currNumOfLevels = prevNumOfLevels + 1;
        end
    end
    
    if prevNumOfLevels ~= currNumOfLevels
        
        startDist = D_Linear(i - 1);
        startENR = 2*energyAlloc(currNumOfLevels)*snrLin(i);
        
    end
    
    if currNumOfLevels == 1
        D_Linear(i) = 1/(1 + 2*energyAlloc(1)*snrLin(i));
    else
        D_Linear(i) = ((1 + startENR)/startENR) * (startDist/(1 + 2*energyAlloc(currNumOfLevels)*snrLin(i)));
%         D_Linear(i) = (startDist/(1 + 2*energyAlloc(currNumOfLevels)*snrLin(i))) * 1 / ((startENR)/(1 + startENR) + (startDist/D_Linear(i-1))*1/(1 + 2*energyAlloc(currNumOfLevels)*snrLin(i)));
    end
    totalEnergy(i) = sum(energyAlloc(1:currNumOfLevels));
    if mod(i,5) == 0
        disp(strcat('Finished SNR = ',num2str(SNR(i))));
    end
end
disp(strcat('Total Energy = ',num2str(sum(energyAlloc(1:currNumOfLevels)))));

if saveResults
    save(strcat('LinearInfDim_Profile',num2str(profileOrder)),'D_Linear','totalEnergy','SNR');
end

figure;
semilogy(SNR,D_Linear,'LineWidth',2.5); hold on;
semilogy(SNR,profile,'-.','LineWidth',2.5);
xlabel('ENR [dB]','FontSize',14); ylabel('Distortion [dB]','FontSize',14);
lgd = legend({'Linear + MLM','Profile'},'FontSize',16,'TextColor','Black');
grid on; grid minor;
title(strcat('Distortion and Profile, order = ',num2str(profileOrder)),'FontSize',14);

