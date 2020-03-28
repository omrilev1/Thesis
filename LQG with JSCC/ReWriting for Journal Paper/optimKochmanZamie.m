clear; close all; clc;

%% Optimize Kochman Zamir JSCC scheme

snr = 7 : 1 : 30;
snrLin = 10.^(snr/10);

Delta = sqrt(12);
rho = [0.6 0.75 0.9 0.99];

MSE = zeros(length(rho),length(snr));
gammaOpt  = zeros(length(rho),length(snr));
gamma = 0.1 : 0.01 : 4;
figure; hold all
lineStyle = ['-g','--k',':c',':g','b',':m'];

for i=1:length(rho)
    
    for j = 1 : length(snr)
        
        MSEanom  = 2*qfunc(0.5 * Delta ./ sqrt(gamma.^2 * (1 - rho(i)^2) + 1/snrLin(j)));
        MSEanom2 = 2*qfunc(3 * 0.5 * Delta ./ sqrt(gamma.^2 * (1 - rho(i)^2) + 1/snrLin(j)));

        MSEweak = (1 - MSEanom) .* ((1 - rho(i)^2)./(1 + gamma.^2 * (1 - rho(i)^2)*snrLin(j)));
        
        MSEtotal = Delta * (gamma*(1 - rho(i)^2) ./ (gamma.^2 * (1 - rho(i)^2) + 1/snrLin(j))) .* MSEanom + ...
            2*Delta * (gamma*(1 - rho(i)^2) ./ (gamma.^2 * (1 - rho(i)^2) + 1/snrLin(j))) .* MSEanom2 + MSEweak;
        [MSEopt,optIdx] = min(MSEtotal);
        gammaOpt = gamma(optIdx);
        MSE(i,j) = MSEopt;
    end
    
    plot(snr,10*log10(1./MSE(i,:)),lineStyle(i),'LineWidth',2)
end
grid on; grid minor;
legend('\rho = 0.6','\rho = 0.75','\rho = 0.9','\rho = 0.99');
xlabel('snr [dB]'); ylabel('SDR [dB]');
title('Optimized MSE vs SNR - Kochman Zamir');