%% Analog PPM Performance simulation
close all; clear all; clc;

ENR = (4:1:12);
ENRlin = 10.^(ENR/10);

Nrun = 2*1e5;
final_MSE = zeros(1,length(ENR));
final_Beta = zeros(1,length(ENR));

for i=1:length(ENR)
    % pick beta to find the optimal point of the curve
    beta_opt = (312*sqrt(pi))^(1/3) * (ENRlin(i))^(-5/6) * exp(ENRlin(i)/6);
    beta = max(1,beta_opt*(0.3125:(1/64):1));
    MSE = zeros(size(beta));% Mean = zeros(size(ENR));
    for beta_idx = 1:length(beta)
        curr_beta = beta(beta_idx);
        W = curr_beta/2;
        
        dt = 1/(100*curr_beta); Fs = 1/dt;
        t = -1.125:dt:1.125;
        tIdx = find(abs(t) <= 0.5);
        
        ppmPulse = zeros(size(t));
        ppmPulse(abs(t) < 1/(2*curr_beta)) = sqrt(curr_beta);
        ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));
        
        Lx = length(ppmPulse);
        Ly = length(ppmPulse)+length(ppmPulse)-1;
        Ly2 = pow2(nextpow2(Ly));
        PPMfreq = fft(ppmPulse, Ly2);
        
        currMSE = zeros(1,Nrun);
        currENR = ENRlin(i);
        parfor n=1:Nrun
            
            % generate source
            S = rand - 0.5;
            
            % PPM Modulation
            TxPulse = zeros(size(t));
            TxPulse(abs(t - S) < 1/(2*curr_beta)) = sqrt(curr_beta);
            TxPulse = sqrt(currENR)*TxPulse/sum(abs(TxPulse.^2)*dt);
            
            noise = randn(size(t));
            noise = sqrt(1/(2*dt))*noise;% sqrt(2*Fs/W)*noise;
            r = TxPulse + noise;
            
            % PPM Correlator receiver
            PPMcorr = fconv(r,Lx,Ly,Ly2,PPMfreq);
            
            %         PPMcorr2 = conv(r,ppmPulse,'same');
            PPMcorr = PPMcorr(tIdx);
            
            [~,maxIdx] = max(PPMcorr);
            sHat = t(min(tIdx) + maxIdx - 1);
            if sHat > 0.5
                sHat = 0.5;
            elseif sHat <= -0.5
                sHat = -0.5;
            end
            
            currMSE(n) = (S - sHat)^2;
        end
        
        MSE(beta_idx) = sum(currMSE)/Nrun;
    end
    [optVal,optIdx] = min(abs(MSE - 13/8 ./ (beta.*ENRlin(i)).^2));
    final_MSE(i) = optVal;
    final_Beta(i) = beta(optIdx);
    
    disp(strcat('Finished ENR = ',num2str(ENR(i))));
end

optSDR_Analytic = 0.058 * (ENRlin).^(-1/3) .* exp(-ENRlin/3);
optBeta_Analytic = (312*sqrt(pi))^(1/3) * (ENRlin).^(-5/6) .* exp(ENRlin/6);

figure;subplot(211);
hold all; 
plot(ENR,10*log10(1./final_MSE),'-.ko','LineWidth',1.5);
plot(ENR,10*log10(1./optSDR_Analytic),'--ms','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empiric Opt SDR','Lower Bound Opt SDR'},'Location','northwest','FontSize',12);
title('Optimal SDR [dB]');
grid on; grid minor;
subplot(212);
semilogy(ENR,final_Beta,'-.ko','LineWidth',1.5);hold on;
semilogy(ENR,optBeta_Analytic,'--ms','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('Optimal Beta [log]');legend({'Empiric Optimal \beta','Optimal \beta - Bound'},'Location','northwest','FontSize',12);
title('Optimal \beta [log]');
grid on; grid minor;



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
