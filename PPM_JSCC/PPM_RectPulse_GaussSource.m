%% Analog PPM Performance simulation
close all; clear ; clc;

ENR = (-1.5:0.5:20);
ENRlin = 10.^(ENR/10);

ENR_opt = 12.5;
ENR_opt_lin = 10^(ENR_opt/10);

c_opt = ((144/pi^2) * sqrt(2*pi*ENR_opt_lin)/(1 + 2*ENR_opt_lin))^(1/3);
beta = c_opt * exp(ENR_opt_lin/6);
W = beta/2;

dt = 1/(500*beta); Fs = 1/dt;

overload = 6.25;
t = -overload:dt:overload;

MSE_MAP = zeros(size(ENR));% Mean = zeros(size(ENR));
MSE_ML = zeros(size(ENR));% Mean = zeros(size(ENR));
 
Nrun =  1e5;

ppmPulse = zeros(size(t));
ppmPulse(abs(t) < 1/(2*beta)) = sqrt(beta);
ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));
% ppmPulse(abs(t) > 0.5) = 0;

Lx = length(ppmPulse);
Ly = length(ppmPulse)+length(ppmPulse)-1;
Ly2 = pow2(nextpow2(Ly));
PPMfreq = fft(ppmPulse, Ly2);

Hist = zeros(length(ENR),Nrun);
TxEnergy = zeros(1,Nrun);
for i=1:length(ENR)
    
    currMSE_MAP = zeros(1,Nrun);currMSE_ML = zeros(1,Nrun);
    currENR = ENRlin(i);
    parfor n=1:Nrun
        
        % generate source - Gaussian Source with edge truncation
        S = randn;
        if abs(S) > max(t)
            S = max(t)*sign(S);
        end
        
        TxPulse = zeros(size(t));
        TxPulse(abs(t - S) < 1/(2*beta)) = sqrt(beta);
        TxPulse = sqrt(currENR)*TxPulse/sum(abs(TxPulse.^2)*dt);
        
        noise = randn(size(t));
        noise = sqrt(1/(2*dt))*noise;
        r = TxPulse + noise;
        
        % PPM Correlator receiver
        PPMcorr = fconv(r,Lx,Ly,Ly2,PPMfreq);

        [~,maxIdx] = max(sqrt(currENR)*PPMcorr*dt - 0.5*(t.^2));
        sHat_MAP = t(maxIdx);
        
        [~,maxIdx] = max(PPMcorr*dt);
        sHat_ML = t(maxIdx);
        
        currMSE_MAP(n) = (S - sHat_MAP)^2;
        currMSE_ML(n) = (S - sHat_ML)^2;
        
    end
    
    MSE_MAP(i) = sum(currMSE_MAP)/Nrun;
    MSE_ML(i) = sum(currMSE_ML)/Nrun;
    
    if mod(i,5) == 0
        disp(strcat('Finished ENR = ',num2str(ENR(i))));
    end
end

figure;hold all;
% plot(ENR,10*log10(min(10,(2*(2*overload)*beta)*(1/(4*sqrt(pi))).*sqrt(ENRlin).*exp(-ENRlin/2).*(1 + sqrt(3/4)*exp(-ENRlin/6)./sqrt(ENRlin) + 4*sqrt(pi./ENRlin)))),'LineWidth',3);
% plot(ENR,10*log10((26/16)./((beta*ENRlin).^(2)) + min(1,(2*(2*overload)*beta)*(1/(4*sqrt(pi))).*sqrt(ENRlin).*exp(-ENRlin/2).*(1 + sqrt(3/4)*exp(-ENRlin/6)./sqrt(ENRlin) + 4*sqrt(pi./ENRlin)))),'LineWidth',3);
plot(ENR,10*log10(min(1.25,(26/16)./((beta*ENRlin).^(2)) + (2*(2*overload)*beta)*(1/(4*sqrt(pi))).*sqrt(ENRlin).*exp(-ENRlin/2).*(1 + sqrt(3/4)*exp(-ENRlin/6)./sqrt(ENRlin) + 4*sqrt(pi./ENRlin)))),'LineWidth',3);
plot(ENR,10*log10(MSE_MAP),'-.','LineWidth',3);
plot(ENR,10*log10(MSE_ML),'--','LineWidth',3);
% plot(ENR,10*log10((26/16)./((beta*ENRlin).^2)),'LineWidth',3);
grid on; grid minor;
xlabel('ENR [dB]'); ylabel('MSE [dB]');
title(strcat('Estimator MSE, ENR_{opt} = ',num2str(ENR_opt)));
legend('Upper Bound','MAP MSE [dB]','ML MSE [dB]');
% ylim([min(10*log10(MSE_MAP))-3 15]);

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
