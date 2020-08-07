%% Analog PPM Performance simulation
close all; clear all; clc;

ENR = 0:0.25:22.5;
ENRlin = 10.^(ENR/10);

ENR_opt = 10;
ENR_opt_lin = 10^(ENR_opt/10);

c_opt = ((144/pi^2) * sqrt(2*pi*ENR_opt_lin)/(2*ENR_opt_lin))^(1/3);
beta = c_opt * exp(ENR_opt_lin/6);
W = beta/2;

dt = 1/(100*beta); Fs = 1/dt;
t = -1.125:dt:1.125;
tIdx = find(abs(t) <= 0.5);

MSE = zeros(size(ENR));% Mean = zeros(size(ENR));

Nrun = 1e4;

ppmPulse = sqrt(2*W)*sin(2*pi*W*t)./(2*pi*W*t);
ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));
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
        
        TxPulse = sqrt(2*W)*sin(2*pi*W*(t - S))./(2*pi*W*(t - S));
        TxPulse(t == S) = max(ppmPulse);
        TxPulse = sqrt(ENRlin(i))*TxPulse/sum(abs(TxPulse.^2)*dt);
        
        noise = randn(size(t));
        noise = sqrt(1/(2*dt))*noise;
        r = TxPulse + noise;
        
        % PPM Correlator receiver
        PPMcorr = fconv(r,Lx,Ly,Ly2,PPMfreq);

%         PPMcorr = conv(r,ppmPulse,'same');

        PPMcorr = PPMcorr(tIdx);
        
        [~,maxIdx] = max(PPMcorr);
        sHat = t(min(tIdx) + maxIdx - 1);
        
        currMSE = currMSE + (S - sHat)^2;
        
%         Hist(i,n) = S - sHat;
    end
    
    MSE(i) = currMSE/Nrun;
    
    if mod(i,5) == 0
        disp(strcat('Finished ENR = ',num2str(ENR(i))));
    end
end

dx = 1e-5;
x = -6.5:dx:6.5;
Pe = zeros(size(ENRlin));
rho = 1; % Exact Term: 4.5
for i=1:length(ENRlin)
    Pe(i) = 1 - (1/sqrt(2*pi))*dx*sum(((1 - qfunc(x + sqrt(2*ENRlin(i)))).^((W-1)*rho)).*exp(-x.^2/2));
end
% analyticTerm = (1./(2*ENRlin))*(3/(pi^2))*(1/beta^2) + ((beta - 1)./sqrt(2*pi*ENRlin)).*exp(-ENRlin/2) .* (1/12);
analyticTerm = (1 - Pe).*(1./(2*ENRlin))*(3/(pi^2))*(1/beta^2) + Pe.* (1/6);

% ZZLB
deps = 1e-4;
eps = 0:deps:1;
eps(eps == 0) = deps;
zzlb = zeros(size(ENRlin));
for i=1:length(ENRlin)
    zzlb(i) = 2 * deps * sum((eps).*(1-eps).*qfunc(sqrt(ENRlin(i) * (1 - sinc(sqrt(2)*beta*eps)))));
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
    J(idx1) = (0.5*(eps1.^2).*((1 - eps1).^2).*exp(-0.5*ENRlin(i)*(1 - sinc(2*beta*eps1))))./(1 - eps1 - (1 - 2*eps1).*exp(-0.25*ENRlin(i)*(1 - sinc(4*beta*eps1))));
    J(idx2) = 0.5*(eps2.^2).*(1 - eps2).*exp(-(ENRlin(i)/2)*(1 - sinc(2*beta*eps2)));
    J = 2*J;
    [wwb(i),minIdx] = max(J);
    optH(i) = eps(minIdx);
end

% % ZZLB Improved
% zzlb_imp = zeros(size(ENRlin));
% dx = 1e-3;
% x = -6.5:dx:6.5;
% deps = 1e-4;
% eps = 0:deps:0.5;
% for i=1:length(ENRlin)
%     
%     % calc curr Pe
%     CurrPe = 1 - (1/sqrt(2*pi))*dx*sum(((1 - qfunc(x + sqrt(2*ENRlin(i)))).^((floor((W)*eps(:))))).*exp(-x.^2/2),2);
%     CurrPe = CurrPe.';
% %     zzlb_imp(i) = 4 * deps * sum((eps).*(1-eps).*max(exp(-ENRlin(i)/2)*(beta*eps - 1)/sqrt(1 + 2*pi*ENRlin(i)),...
% %         qfunc(sqrt(ENRlin(i) * (1 - sinc(2*beta*eps))))));
% 
%     zzlb_imp(i) = 4 * deps * sum((eps).*(1-eps).*max(0.5*CurrPe,qfunc(sqrt(ENRlin(i) * (1 - sinc(2*beta*eps))))));
% end

% Analytic upper bound
deps = 1e-5;
eps = 0:deps:1;
eps(eps == 0) = deps;
upperBound = zeros(size(ENRlin));
for i=1:length(ENRlin)
    currNorm = deps*sum(qfunc(sqrt(ENRlin(i) * (1 - sinc(2*beta*eps)))));
    upperBound(i) = 2 * deps * sum((eps.^2).*(1-eps).*qfunc(sqrt(ENRlin(i) * (1 - sinc(2*beta*eps))))/currNorm);
end
figure;hold all;
plot(ENR,10*log10(upperBound),'--','LineWidth',3);
% plot(ENR,10*log10(analytic2_imp),'LineWidth',3);
plot(ENR,10*log10(MSE),'-*','LineWidth',3);
plot(ENR,10*log10(analyticTerm),'-.','LineWidth',3);
% plot(ENR,10*log10((1/6)./(1 + ENRlin)),'LineWidth',3);
plot(ENR,10*log10(wwb),'LineWidth',3);
plot(ENR,10*log10(zzlb),'LineWidth',3);

grid on; grid minor;
xlabel('ENR [dB]'); ylabel('MSE [dB]');
title(strcat('Estimator MSE, ENR_{opt} = ',num2str(ENR_opt))); 
legend('Upper Bound','Empiric','Wozencraft-Jacobs','WWB','ZZLB');

% legend('Upper Bound','Empiric','Wozencraft-Jacobs','zzlb imp','zzlb');
% legend('Upper Bound','Upper Bound Improved','Empiric','Wozencraft-Jacobs','Linear Transmission');


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