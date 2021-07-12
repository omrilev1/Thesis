%% Simulation of Successive MLM and Linear transmissions, for the transmissions of an uniform source
%% over a channel with unlimited bandwidth. We optimize the scheme parameters according to the simulated profile

close all; clear all; clc;

profileOrder = 2;
saveResults = 1;
 
% Init parameters and arrays structures
SNR = -7:0.25:35;  % 8              % Direct channel SNR
snrLin = 10.^(SNR/10);

profile = 1./(1 + (snrLin).^profileOrder);

maxStages = 10; % number of PPM levels
beta = zeros(1,maxStages);

%% Modulo parameters:
% interval
IntLin = sqrt(12)*ones(1,maxStages); % [sqrt(12)/2.125 sqrt(12)/2.125 sqrt(12)/2.125];

Delta = 0.9; alpha = 0.675;
energyAlloc = ones(1,maxStages) .* Delta .* exp(-(0:1:(maxStages - 1)) * alpha);

%% Simulation parameters
N_avg = 2^13;
E_0 = 1;

D_PPM = zeros(length(SNR),1);
totalEnergy = zeros(length(SNR),1);

currNumOfLevels = 1;finished = 0;
% currE = zeros(1,currNumOfLevels);currE(1) = 1;

etaPPM = 1;

for i=1:length(SNR)
    
    prevNumOfLevels = currNumOfLevels;
    if i > 1
        if D_PPM(i - 1) > profile(i - 1)
            currNumOfLevels = prevNumOfLevels + 1;
            
        end
    end
    
    if prevNumOfLevels ~= currNumOfLevels
        
        % calculate optimal ppm beta and generate ppm pulse
        ENR_opt_lin = 10^(-2.5/10) * energyAlloc(currNumOfLevels) * 10^(SNR(i)/10);
        
        beta(currNumOfLevels) = 0.5 * (312 * sqrt(pi)) ^ (1/3) * exp(ENR_opt_lin/6) .* (ENR_opt_lin .^(-5/6));
        
        dt = 1/(650*beta(currNumOfLevels));
        t = -1.125:dt:1.125;
        tIdx = find(abs(t) <= 0.5);
        
        ppmPulse = zeros(size(t));
        ppmPulse(abs(t) < 1/(2*beta(currNumOfLevels))) = sqrt(beta(currNumOfLevels));
        ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));
        
        Lx = length(ppmPulse);
        Ly = length(ppmPulse)+length(ppmPulse)-1;
        Ly2 = pow2(nextpow2(Ly));
        PPMfreq{currNumOfLevels - 1} = fft(ppmPulse, Ly2);
        
        D_PPM_UB = 0.058 * (ENR_opt_lin).^(-1/3) .* exp(-ENR_opt_lin/3);
        
        currEtaPPM = optimizeEta(D_PPM(i-1),IntLin(currNumOfLevels),D_PPM_UB);
        currEtaPPM = 1.35 * currEtaPPM;
        etaPPM = [etaPPM currEtaPPM];
        
    end
    
    currDist_PPM = zeros(1,N_avg);
    currE_tot = zeros(1,N_avg);
    for n = 1 : N_avg
        
        % calculate the sources of the current transmission levels
        S = rand - 0.5;
        TxSourcePpm = zeros(1,maxStages);
        dither = zeros(1,maxStages);
        
        TxSourcePpm(1) = S;
        
        if currNumOfLevels > 1
            for idx=2:currNumOfLevels
                delta = IntLin(idx);
                TxSourcePpm(idx) = mod(etaPPM(idx)*S + delta,2*delta) - delta;
            end
        end
        
        % simulate the successive transmission and reception
        y_PPM = zeros(1,currNumOfLevels);
        
        sHat_PPM = 0;
        for k=1:currNumOfLevels
            
            %% channel part
            if k==1
                
                y_PPM(1) = sqrt(12) * sqrt(energyAlloc(k)) * S + sqrt(1/snrLin(i))*randn;
                
            else
                
                dt = 1/(650*beta(k));
                t = -1.125:dt:1.125;
                tIdx = find(abs(t) <= 0.5);
                
                Lx = length(t);
                Ly = length(t)+length(t)-1;
                Ly2 = pow2(nextpow2(Ly));
                
                TxPulse = zeros(size(t));
                TxPulse(abs(t - TxSourcePpm(k) / (2*delta)) < 1/(2*beta(k))) = sqrt(beta(k));
                TxPulse = sqrt(energyAlloc(k)*snrLin(i)) * TxPulse/sum(abs(TxPulse.^2)*dt);
                
                noise = randn(size(t));
                noise = sqrt(1/(2*dt))*noise;% sqrt(2*Fs/W)*noise;
                r = TxPulse + noise;
                
                % PPM Correlator receiver
                PPMcorr = fconv(r,Lx,Ly,Ly2,PPMfreq{k - 1});
                PPMcorr = PPMcorr(tIdx);
                
                [~,maxIdx] = max(PPMcorr);
                y_PPM(k) = t(min(tIdx) + maxIdx - 1);
                if y_PPM(k) > 0.5
                    y_PPM(k) = 0.5;
                elseif y_PPM(k) <= -0.5
                    y_PPM(k) = -0.5;
                end
                y_PPM(k) = y_PPM(k) * 2 * delta;
            end
            
            %% MLM Part
            if k==1
                sHat_PPM = y_PPM(1) / (sqrt(12) * sqrt(energyAlloc(k))) ;
            else
                sHat_PPM = sHat_PPM + (1/etaPPM(k)) * (mod(y_PPM(k) - etaPPM(k)*sHat_PPM + IntLin(k-1),2*IntLin(k-1)) - IntLin(k-1));
            end
            
        end
        
        currDist_PPM(n) = (sHat_PPM - S)^2;
        currE_tot(n) = sum(energyAlloc);
    end
    
    D_PPM(i) = sum(currDist_PPM)/N_avg;
    totalEnergy(i) = sum(energyAlloc(1:currNumOfLevels));
    
    if mod(i,5) == 0
        disp(strcat('Finished SNR = ',num2str(SNR(i))));
    end
    
end

if saveResults
    save(strcat('PPM_Profile',num2str(profileOrder)),'D_PPM','totalEnergy','Delta','alpha','SNR','profile');
end

figure;
semilogy(SNR,D_PPM,'LineWidth',2.5); hold on;
semilogy(SNR,profile,'-.','LineWidth',2.5);
xlabel('ENR [dB]','FontSize',14); ylabel('Distortion [dB]','FontSize',14);
lgd = legend({'PPM + MLM','Profile'},'FontSize',16,'TextColor','Black');
grid on; grid minor;
title(strcat('Distortion and Profile, order = ',num2str(profileOrder)),'FontSize',14);


function etaOpt = optimizeEta(D,Int,sigmaN)
eta = 1:0.0125:200000;

% calculate probability of error
P_e = 2 * qfunc(Int./sqrt(eta.^2*D + sigmaN));

% calculate distortion
dist = (1 - P_e).*(sigmaN./(eta.^2)) + P_e.*(2*Int^2);

[~,etaOpt_idx] = min(dist);
etaOpt = eta(etaOpt_idx);
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

