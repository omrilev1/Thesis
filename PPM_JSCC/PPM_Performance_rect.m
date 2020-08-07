%% Analog PPM Performance simulation
close all; clear all; clc;

ENR = (-3:0.25:18);
ENRlin = 10.^(ENR/10);

ENR_opt = 10;
ENR_opt_lin = 10^(ENR_opt/10);

c_opt = ((144/pi^2) * sqrt(2*pi*ENR_opt_lin)/(1 + 2*ENR_opt_lin))^(1/3);
beta = c_opt * exp(ENR_opt_lin/6);
W = beta/2;

dt = 1/(100*beta); Fs = 1/dt;
t = -1.125:dt:1.125;
tIdx = find(abs(t) <= 0.5);

MSE = zeros(size(ENR));% Mean = zeros(size(ENR));

Nrun = 1e4;

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
    currMSE = 0;
    %     currMean = 0;
    
    for n=1:Nrun
        
        % generate source
        S = rand - 0.5;
         
        TxPulse = zeros(size(t));
        TxPulse(abs(t - S) < 1/(2*beta)) = sqrt(beta);
        TxPulse = sqrt(ENRlin(i))*TxPulse/sum(abs(TxPulse.^2)*dt);
        
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
        
        currMSE = currMSE + (S - sHat)^2;
        Hist(i,n) = S - sHat;
    end
    
    MSE(i) = currMSE/Nrun;
    
    if mod(i,5) == 0
        disp(strcat('Finished ENR = ',num2str(ENR(i))));
    end
end

% analytic 2
deps = 1e-6;
eps = 0:deps:1;
upperBound = zeros(size(ENRlin));
autoCorr = zeros(size(eps));
autoCorr(eps < 1/(beta)) = 1 - beta*eps(eps < 1/(beta));
for i=1:length(ENRlin)
        currNorm = deps*sum(qfunc(sqrt(ENRlin(i) * (1 - autoCorr))));
    upperBound(i) = 2 * deps * sum((eps.^2).*(1-eps).*qfunc(sqrt(ENRlin(i) * (1 - autoCorr)))/currNorm);
    
end

% ZZLB
deps = 1e-6;
eps = 0:deps:1;
eps(eps == 0) = deps;
zzlb = zeros(size(ENRlin));
autoCorr = zeros(size(eps));
autoCorr(eps < 1/(beta)) = 1 - beta*eps(eps < 1/(beta));
for i=1:length(ENRlin)
    zzlb(i) = 2 * deps * sum((eps).*(1-eps).*qfunc(sqrt(ENRlin(i) * (1 - autoCorr))));
end

% WWB
deps = 1e-5;
eps = 0:deps:1;
eps(eps == 0) = deps;
wwb = zeros(size(ENRlin));optH = zeros(size(ENRlin));
for i=1:length(ENRlin)
    J = zeros(size(eps));
    idx1 = (eps <= 1/2); idx2 = (eps > 1/2);
    eps1 = eps(idx1); eps2 = eps(idx2);
    autoCorr = zeros(size(eps1)); autoCorr(eps1 < 1/(beta)) = 1 - beta*eps1(eps1 < 1/(beta));
    autoCorr2 = zeros(size(eps1)); autoCorr2(eps1 < 1/(beta/2)) = 1 - (beta/2)*eps1(eps1 < 1/(beta/2));
    
    J(idx1) = (0.25*(eps1.^2).*((1 - eps1).^2).*exp(-0.5*ENRlin(i)*(1 - autoCorr)))./(1 - eps1 - (1 - 2*eps1).*exp(-0.25*ENRlin(i)*(1 - autoCorr2)));
    J(idx2) = 0.5*(eps2.^2).*(1 - eps2).*exp(-(ENRlin(i)/2));

    [wwb(i),minIdx] = max(J);
    optH(i) = eps(minIdx);
end

figure;hold all;
plot(ENR,10*log10((26/16)./((beta*ENRlin).^2) + min(1/6,(beta/24)*(1/(4*sqrt(pi))).*sqrt(ENRlin).*exp(-ENRlin/2).*(1 + sqrt(3/4)*exp(-ENRlin/6)./sqrt(ENRlin) + 4*sqrt(pi./ENRlin)))),':','LineWidth',3);
% plot(ENR,10*log10(min(1/6,(beta/24)*(1/(4*sqrt(pi))).*sqrt(ENRlin).*exp(-ENRlin/2).*(1 + sqrt(3/4)*exp(-ENRlin/6)./sqrt(ENRlin) + 4*sqrt(pi./ENRlin)))),'LineWidth',3);
plot(ENR,10*log10(MSE),'-.','LineWidth',3);
plot(ENR,10*log10(zzlb),'.','LineWidth',3);
plot(ENR,10*log10(wwb),':','LineWidth',3);
% plot(ENR,10*log10((26/16)./((beta*ENRlin).^2)),'LineWidth',3);
% plot(ENR,10*log10((1/6)*exp(-ENRlin/2)),'LineWidth',3);

grid on; grid minor;
xlabel('ENR [dB]'); ylabel('MSE [dB]');
title(strcat('Estimator MSE, ENR_{opt} = ',num2str(ENR_opt)));
legend('Upper Bound','Empiric','ZLLB','WWB');
% legend('Upper Bound','Empiric','Linear Transmission');
ylim([min(10*log10(zzlb)) , 0])
xlim([-5,18])

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
