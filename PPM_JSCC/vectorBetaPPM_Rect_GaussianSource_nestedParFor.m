%% Analog PPM Performance simulation
close all; clear all; clc;

ENR = (4:1:15);
ENRlin = 10.^(ENR/10);

Nrun = 1e5;
final_MSE = zeros(1,length(ENR));
final_Beta = zeros(1,length(ENR));

for i=1:length(ENR)
    % pick beta to find the optimal point of the curve
    beta_opt = (13/8)^(1/3) * (ENRlin(i))^(-5/6) * exp(ENRlin(i)/6);
    beta = beta_opt*(0.5:(1/8):3.5);
    MSE = zeros(size(beta));
    currENRlin = ENRlin(i);
    parfor beta_idx = 1:length(beta)
        curr_beta = beta(beta_idx);
        currMSE = simulateGaussianPPM(curr_beta,currENRlin,Nrun,6)
        MSE(beta_idx) = currMSE;
    end
    %     [optVal,optIdx] = min(abs(MSE - 13/8 ./ (beta.*currENRlin).^2));
    [optVal,optIdx] = min(MSE);
    final_MSE(i) = MSE(optIdx);
    final_Beta(i) = beta(optIdx);
    disp(strcat('Finished ENR = ',num2str(ENR(i))));
end
optSDR_Analytic_Approx = 3 * (13/8)^(1/3) * (ENRlin).^(-1/3) .* exp(-ENRlin/3);
optBeta_Analytic = (13/8)^(1/3) .* (ENRlin).^(-5/6) .* exp(ENRlin/6);

% exact upper bound

D_S = ((13/8) + 2 * (sqrt(2*optBeta_Analytic.*ENRlin) - 1).*exp(-ENRlin .* (1 - 1./sqrt(2*optBeta_Analytic.*ENRlin)).^2)) ./ ((sqrt(optBeta_Analytic.*ENRlin) - 1/sqrt(2)).^4) ...
    + exp(-optBeta_Analytic.*ENRlin)./(optBeta_Analytic.^2);
D_L = 2*optBeta_Analytic.*sqrt(ENRlin).*exp(-ENRlin/2) .* (1 + 3*sqrt(2*pi./ENRlin) + 12*exp(-1)./(optBeta_Analytic.*sqrt(ENRlin)) ...
    + 8*exp(-1)./(sqrt(8*pi)*optBeta_Analytic) + sqrt(8./(pi*ENRlin)) + 12^(3/2) * exp(-3/2) ./(optBeta_Analytic.*sqrt(32*pi*ENRlin)));
optSDR_Analytic_Approx = (D_S + D_L);
optSDR_Analytic_Approx = min(optSDR_Analytic_Approx,1);

% optimize the lower bound
% exact upper bound
optSDR_Analytic_Exact = zeros(size(ENR));
opt_beta = zeros(size(ENR));
for i=1:length(ENR)
    beta_opt = (13/8)^(1/3) .* (ENRlin(i)).^(-5/6) .* exp(ENRlin(i)/6);
    beta_vec =  max(1,beta_opt*(0.0125:(1/128):5));
    
    D_t = ((13/8) + 2 * (sqrt(2*beta_vec.*ENRlin(i)) - 1).*exp(-ENRlin(i) .* (1 - 1./sqrt(2*beta_vec.*ENRlin(i))).^2)) ./ ((sqrt(beta_vec.*ENRlin(i)) - 1/sqrt(2)).^4) ...
        + exp(-beta_vec.*ENRlin(i))./(beta_vec.^2);
    
    D_L = 2*beta_vec.*sqrt(ENRlin(i)).*exp(-ENRlin(i)/2) .* (1 + 3*sqrt(2*pi./ENRlin(i)) + 12*exp(-1)./(beta_vec.*sqrt(ENRlin(i))) ...
        + 8*exp(-1)./(sqrt(8*pi)*beta_vec) + sqrt(8./(pi*ENRlin(i))) + 12^(3/2) * exp(-3/2) ./(beta_vec.*sqrt(32*pi*ENRlin(i))));
    
    [optSDR_Analytic_Exact(i),optIdx] = min(D_t + D_L);
    opt_beta(i) = beta_vec(optIdx);
end

figure;subplot(211);
hold all;
plot(ENR,10*log10(1./final_MSE),'-.ko','LineWidth',1.5);
plot(ENR,10*log10(1./optSDR_Analytic_Exact),'--ms','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empiric Opt SDR','Lower Bound Opt SDR'},'Location','northwest','FontSize',12);
grid on; grid minor;
subplot(212);
semilogy(ENR,final_Beta,'-.ko','LineWidth',1.5);hold on;
semilogy(ENR,opt_beta,'--ms','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('Optimal Beta [log]');legend({'Empiric Optimal \beta','Optimal \beta - Bound'},'Location','northwest','FontSize',12);
grid on; grid minor;

figure;
hold all;
plot(ENR,10*log10(1./final_MSE),'-.ko','LineWidth',1.5);
plot(ENR,10*log10(1./optSDR_Analytic_Approx),'--ms','LineWidth',1.5);
plot(ENR,10*log10(1./optSDR_Analytic_Exact),'->','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empiric Opt SDR','Lower Bound Opt SDR - Approx','Lower Bound Opt SDR - Exact'},'Location','northwest','FontSize',12);
title('Optimal SDR [dB]');
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

function [MSE] = simulateGaussianPPM(beta,SNRlin,Nrun,overload)

dt = 1/(250*beta);

t = -overload:dt:overload;

ppmPulse = zeros(size(t));
ppmPulse(abs(t) < 1/(2*beta)) = sqrt(beta);
ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));
% ppmPulse(abs(t) > 0.5) = 0;

Lx = length(ppmPulse);
Ly = length(ppmPulse)+length(ppmPulse)-1;
Ly2 = pow2(nextpow2(Ly));
PPMfreq = fft(ppmPulse, Ly2);
currMSE_MAP = zeros(1,Nrun);

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

end

MSE = sum(currMSE_MAP)/Nrun;


end

