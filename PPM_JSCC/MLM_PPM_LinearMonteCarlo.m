%% Simulation of Successive MLM and Linear transmissions, for the transmissions of an uniform source 
%% over a channel with unlimited bandwidth. We optimize the scheme parameters according to the simulated profile 

close all; clear all; clc;

profileOrder = 2;
saveResults = 0;

% Init parameters and arrays structures
SNR = -7:0.25:35;  % 8              % Direct channel SNR
snrLin = 10.^(SNR/10);

profile = 1./(1 + (snrLin).^profileOrder);

maxStages = 40; % number of PPM levels

%% Modulo parameters:
% interval
IntLin = sqrt(12)*ones(1,maxStages); % [sqrt(12)/2.125 sqrt(12)/2.125 sqrt(12)/2.125];

Delta = 0.9; alpha = 0.64;
energyAlloc = ones(1,maxStages) .* Delta .* exp(-(0:1:(maxStages - 1)) * alpha);

%% Simulation parameters
N_avg = 2^14;

D_Linear = zeros(length(SNR),1);
totalEnergy = zeros(length(SNR),1);

currNumOfLevels = 1;finished = 0;
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
        
        currEtaLinear = optimizeEta(D_Linear(i-1),IntLin(currNumOfLevels),1/(1 + energyAlloc(currNumOfLevels) * snrLin(i)));
%         currEnergy = energyAlloc(currNumOfLevels);% calcModuloEnergy(IntLin(currNumOfLevels-1),currEtaLinear);
        etaLinear = [etaLinear currEtaLinear];
%         currE = [currE currEnergy];
        
    end
    currDist_Linear = zeros(1,N_avg);
    currE_tot = zeros(1,N_avg);
    parfor n = 1 : N_avg
        
        % calculate the sources of the current transmission levels
        S = rand - 0.5;
        TxSourceLinear = zeros(1,maxStages);
        dither = zeros(1,maxStages);
        
        TxSourceLinear(1) = S;
        
        if currNumOfLevels > 1
            for idx=2:currNumOfLevels
                delta = IntLin(idx);
                TxSourceLinear(idx) = mod(etaLinear(idx)*S + delta,2*delta) - delta;
            end
            
        end
                
        % simulate the successive transmission and reception
        y_linear = zeros(1,currNumOfLevels);
        
        sHat_Linear = 0;
        sHat_PPM = 0;
        for k=1:currNumOfLevels
            
            %% channel part
            if k==1
                y_linear(1) = sqrt(12) * sqrt(energyAlloc(k)) * S + sqrt(1/snrLin(i))*randn;
            else
                y_linear(k) = sqrt(energyAlloc(k)) * TxSourceLinear(k) + sqrt(1/snrLin(i))*randn;
            end
            
            %% MLM Part
            if k==1
                sHat_Linear = y_linear(1) / (sqrt(12) * sqrt(energyAlloc(k))) ;
            else
                alphaMMSE_Linear = sqrt(energyAlloc(k))/(energyAlloc(k) + 1/snrLin(i));
                sHat_Linear = sHat_Linear + (1/etaLinear(k)) * (mod(alphaMMSE_Linear * y_linear(k) - etaLinear(k)*sHat_Linear + IntLin(k-1),2*IntLin(k-1)) - IntLin(k-1));
            end
                        
        end
        
        currDist_Linear(n) = (sHat_Linear - S)^2;
        currE_tot(n) = sum(energyAlloc);
    end
    totalEnergy(i) = sum(energyAlloc(1:currNumOfLevels));
    D_Linear(i) = sum(currDist_Linear)/N_avg;
   
    if mod(i,5) == 0
        disp(strcat('Finished SNR = ',num2str(SNR(i))));
    end
    
end

if saveResults
    save(strcat('Linear_Profile',num2str(profileOrder)),'D_Linear','totalEnergy','Delta','alpha','SNR');
end

figure;
semilogy(SNR,D_Linear,'LineWidth',2.5); hold on;
semilogy(SNR,profile,'-.','LineWidth',2.5); 
xlabel('ENR [dB]','FontSize',14); ylabel('Distortion [dB]','FontSize',14);
lgd = legend({'Linear + MLM','Profile'},'FontSize',16,'TextColor','Black');
grid on; grid minor;
title(strcat('Distortion and Profile, order = ',num2str(profileOrder)),'FontSize',14);



function etaOpt = optimizeEta(D,Int,sigmaN)
eta = 1:0.0125:300000;

% calculate probability of error
P_e = 2 * qfunc(Int./sqrt(eta.^2*D + sigmaN));

% calculate distortion
dist = (1 - P_e).*(sigmaN./(eta.^2)) + P_e.*(2*Int^2);

[~,etaOpt_idx] = min(dist);
etaOpt = eta(etaOpt_idx);
end
