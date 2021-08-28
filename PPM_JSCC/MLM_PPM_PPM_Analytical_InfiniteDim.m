%% Simulation of Successive MLM and Linear transmissions, for the transmissions of an uniform source
%% over a channel with unlimited bandwidth. We optimize the scheme parameters according to the simulated profile

close all; clear all; clc;

profileOrder = 2;
saveResults = 1;

% Init parameters and arrays structures
SNR = -7:0.25:30;  % 8              % Direct channel SNR
snrLin = 10.^(SNR/10);
profile = 1./(1 + (snrLin).^profileOrder);

% profile = 1./(1 + (snrLin/snrLin(1)).^profileOrder);

maxStages = 20; % number of PPM levels

%% Modulo parameters:

% Delta = 0.85; alpha = 0.55;
% M = 6;
Delta = 0.848; alpha = 0.55;
M = 6;
energyAlloc = ones(1,M + 1) .* Delta .* exp(-(0:1:(M)) * alpha);

disp(strcat('Total Energy = ',num2str(sum(energyAlloc))))
%% Simulation parameters

D_PPM_Anal = zeros(length(SNR),1);
SDR_tot = ones(length(SNR),1);
totalEnergy = ones(length(SNR),1);
currNumOfLevels = 1;
% currE = zeros(1,currNumOfLevels);currE(1) = 1;

etaLinear = 1;

for i=1:length(SNR)
    
    prevNumOfLevels = currNumOfLevels;
    if ((currNumOfLevels <= M))
        if i>4
            if ((D_PPM_Anal(i - 1) > profile(i - 1)))
                currNumOfLevels = prevNumOfLevels + 1;
            end
        end
    end
    
    if prevNumOfLevels ~= currNumOfLevels
        startDist = D_PPM_Anal(i - 1);
        startENR = 2*energyAlloc(currNumOfLevels)*snrLin(i);
        
        if currNumOfLevels > M
            % calculate optimal ppm beta and generate ppm pulse
            beta = 0.9;
            EN0_Start = energyAlloc(currNumOfLevels)*snrLin(i);
            D_S = ((13/8) + sqrt(2/beta) * (sqrt(2*beta*EN0_Start) - 1)*exp(-EN0_Start * (1 - 1/sqrt(2*beta*EN0_Start))^2)) / ((sqrt(beta*EN0_Start) - 1/sqrt(2))^4) ...
                + exp(-beta*EN0_Start)/beta^2;
            D_L = 2*beta*sqrt(EN0_Start)*exp(-EN0_Start/2) * (1 + 3*sqrt(2*pi/EN0_Start) + 12*exp(-1)/(beta*sqrt(EN0_Start)) ...
                + 8*exp(-1)/(sqrt(8*pi)*beta) + sqrt(8/(pi*EN0_Start)) + 12^(3/2) * exp(-3/2) /(beta*sqrt(32*pi*EN0_Start)));
            startSDR = 1 /(D_S + D_L);
            startSDR = max(startSDR,1);
        end
    end
    
    if currNumOfLevels == 1
        D_PPM_Anal(i) = 1/(1 + 2*energyAlloc(1)*snrLin(i));
    else
        
        if currNumOfLevels > M
            % substitute PPM SDR Upper bound for current ENR
            EN0 = energyAlloc(2)*snrLin(i);
            D_S = ((13/8) + sqrt(2/beta) * (sqrt(2*beta*EN0) - 1)*exp(-EN0 * (1 - 1/sqrt(2*beta*EN0))^2)) / ((sqrt(beta*EN0) - 1/sqrt(2))^4) ...
                + exp(-beta*EN0)/beta^2;
            D_L = 2*beta*sqrt(EN0)*exp(-EN0/2) * (1 + 3*sqrt(2*pi/EN0) + 12*exp(-1)/(beta*sqrt(EN0)) ...
                + 8*exp(-1)/(sqrt(8*pi)*beta) + sqrt(8/(pi*EN0)) + 12^(3/2) * exp(-3/2) /(beta*sqrt(32*pi*EN0)));
            SDR_PPM = 1 /(D_S + D_L);
            if SDR_PPM < 1
                SDR_PPM = 1;
            end
            %             D_PPM_Anal(i) = (startDist/(1 + SDR_PPM)) ...
            %                 * 1/(startSDR/(1 + startSDR) + (startDist/D_PPM_Anal(i-1)) * 1/(1 + SDR_PPM));
            
            %             D_PPM_Anal(i) = (startDist/SDR_PPM) ...
            %                         * 1/(1 + (startDist/D_PPM_Anal(i-1)) * (1/SDR_PPM));
            
            D_PPM_Anal(i) = (startDist/SDR_PPM);
            
            SDR_tot(i) = SDR_PPM;
        else
            
            D_PPM_Anal(i) = ((1 + startENR)/startENR) * (startDist/(1 + 2*energyAlloc(currNumOfLevels)*snrLin(i)));
        end
    end
    totalEnergy(i) = sum(energyAlloc(1:currNumOfLevels));
    
    if mod(i,5) == 0
        disp(strcat('Finished SNR = ',num2str(SNR(i))));
    end
end

if saveResults
    save(strcat('PpmInfDim_Anal_Profile',num2str(profileOrder)),'D_PPM_Anal','SNR','totalEnergy');
end

figure;
semilogy(SNR,D_PPM_Anal,'LineWidth',2.5); hold on;
semilogy(SNR,profile,'-.','LineWidth',2.5);
xlabel('ENR [dB]','FontSize',14); ylabel('Distortion [dB]','FontSize',14);
lgd = legend({'PPM + MLM','Profile'},'FontSize',16,'TextColor','Black');
grid on; grid minor;
title(strcat('Distortion and Profile, order = ',num2str(profileOrder)),'FontSize',14);


bypassAnalytical = 1;
%% Analytical calculation:

if ~bypassAnalytical
    Delta = 0.01:0.01:10;
    alpha = 0.01:0.01:10;
    [DeltaGrid,alphaGrid] = meshgrid(Delta,alpha);
    M = 1:1:5;
    optVal = zeros(size(M));
    beta = 0.6;
    
    for i=1:length(M)
        
        % Energy of first M - 1 layers
        if M(i) == 1
            E1 = exp(alphaGrid)./Delta;
        else
            E1 = exp(alphaGrid)./Delta + 0.5*Delta.*(exp(2*alphaGrid) - 1).*(1 + sqrt(1 + 4*exp(3*alphaGrid)./((1 - exp(2*alphaGrid)).^2))).*...
                exp(-2*alphaGrid).*(1 - exp(-alphaGrid*(M(i) - 1)))./(1 - exp(-1*alphaGrid));
        end
        
        % find the energy of the last layer
        N_M = DeltaGrid.*exp(-alphaGrid * M(i));
        E_M = zeros(size(N_M));
        for idx1 = 1:size(N_M,1)
            for idx2 = 1:size(N_M,2)
                Nvec = ((1/100):(1/100):100)*N_M(idx1,idx2);
                profile = (1./(1 + (1./Nvec).^2));
                EN0 = 1./Nvec;
                
                D_S = ((13/8) + 2 * (sqrt(2*beta*EN0) - 1).*exp(-EN0 .* (1 - 1./sqrt(2*beta*EN0)).^2)) ./ ((sqrt(beta*EN0) - 1/sqrt(2)).^4) ...
                    + exp(-beta*EN0)/beta^2;
                D_L = 2*beta*sqrt(EN0).*exp(-EN0/2) .* (1 + 3*sqrt(2*pi./EN0) + 12*exp(-1)./(beta*sqrt(EN0)) ...
                    + 8*exp(-1)/(sqrt(8*pi)*beta) + sqrt(8./(pi*EN0)) + 12^(3/2) * exp(-3/2)./(beta*sqrt(32*pi*EN0)));
                SDR_PPM = 1 ./(D_S + D_L);
                
                [minVal,minIdx] = min(profile - ((1./(1 + (1./N_M(idx1,idx2)).^2)))./SDR_PPM);
                
                if minVal > 0
                    E_M(idx1,idx2) = 1/Nvec(minIdx);
                else
                    E_M(idx1,idx2) = 1e4;
                end
            end
        end
        optVal(i) = min(min(E1 + E_M));
    end
end