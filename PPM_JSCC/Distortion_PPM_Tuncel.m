% close all;
% clear all; clc;

% Infinite BW JSCC schemes

ENR_min = 12.5;
ENR_min_lin = 10^(ENR_min/10);
ENR = -3:0.25:22;
ENR_lin = 10.^(ENR/10);

% Analog PPM
c_opt = ((144/pi^2) * sqrt(2*pi*ENR_min_lin)/(1 + 2*ENR_min_lin))^(1/3);
beta = c_opt * exp(ENR_min_lin / 6);
MSE_Anlog = 12./(pi^2 * beta^2 * (1 + 2*ENR_lin)) + (1/6) * ((beta - 1)./sqrt(2*pi*ENR_lin)) .* exp(-ENR_lin / 2);

% Tuncel optimal + compander
beta = 0.1385; c = 0.8281;
dx = 1e-5;
x = -0.5:dx:0.5;
lambda = 1./(6^(1/3) * c^(2/3) * (2*c*(x.^2) + beta).^(1/3));
lambda = lambda/sum(lambda*dx);
int1 = sum(lambda.*(x.^2)*dx);
int2 = sum(dx./(lambda.^2));
MSE_Tuncel = exp(-ENR_lin/3)*(c/6 + 2*c*int1) + int2 * exp(-ENR_min_lin/3)/(12*c^2);

% ZZLB
% dh = 1e-4;
% h = dh:dh:1;
% beta = c_opt * exp(ENR_min_lin / 6);
% ZZLB = beta * dh*sum(h.*(1-h).*qfunc(sqrt(ENR_lin(:).*(1 - sin(2*pi*beta*h)./(2*pi*beta*h)))),2);
dh = 0.5*1e-6;
h = dh:dh:0.5;
beta = c_opt * exp(ENR_min_lin / 6);
ZZLB = 7.9*dh*sum(h.*(0.5-h).*qfunc(sqrt(ENR_lin(:).*(1 - sin(pi*beta*h)./(pi*beta*h)))),2);

figure; hold all
plot(ENR,10*log10(MSE_Tuncel),'-','LineWidth',2)
plot(ENR,10*log10(MSE_Anlog),'--','LineWidth',2)
plot(ENR,10*log10(ZZLB),'-.','LineWidth',2)
plot(ENR_min*ones(size(ENR)),linspace(min(10*log10(ZZLB)),max(10*log10(MSE_Tuncel)),length(ENR)),'--','LineWidth',2)
grid on; grid minor;
xlabel('\gamma [dB]');
ylabel('MSE [dB]');
title('MSE vs ENR');
legend('Tuncel Separation Scheme + Companding','Analog PPM','Ziv-Zakai Lower Bound','Optimization ENR');


%% Calculate number of levels
figure;hold all;
n_0 = 1:1:4;
ENR_min = 10; ENR_min_lin = 10^(ENR_min/10);
ENR = 5:0.25:60; ENR_lin = 10.^(ENR/10);
c_opt = ((72/pi^2) * sqrt(2*pi*ENR_min_lin)/ENR_min_lin)^(1/3);
beta = c_opt * exp(ENR_min_lin / 6);

plot(ENR,10*log10(1./(1 + 2*ENR_lin.^2)),'LineWidth',2); % desired profile
plot(ENR,10*log10((1/12)*((beta - 1)*exp(-ENR_lin/2))./sqrt(2*pi*ENR_lin) + 6./(beta^2*pi^2*ENR_lin)),'LineWidth',2);

currBound = zeros(length(n_0),length(ENR));
for i=1:length(n_0)
    
    if i==1
        curr_ENR = ENR_lin/5; % ENR_lin/((i+1)^3);
    elseif i==2
        curr_ENR = ENR_lin/75;
    elseif i==3
        curr_ENR = ENR_lin/625;
    elseif i==4
        curr_ENR = ENR_lin/5000;
    end
    currBound(i,:) = (beta^(-2*(i-1))) * exp(-curr_ENR/2) .* (beta - 1)./(12*sqrt(2*pi*curr_ENR)) + beta^(-2*i) * (6./(pi^2*beta^2*curr_ENR) + (1/12)*exp(-curr_ENR/2).*(beta - 1)./sqrt(2*pi*curr_ENR));
    plot(ENR,10*log10(currBound(i,:)),'LineWidth',2);
    
end

grid on; grid minor;
xlabel('ENR [dB]');
ylabel('Distortion [dB]');
legend('D = 1./(1 + \alpha \gamma ^2)','n_0 = 0','n_0 = 1','n_0 = 2','n_0 = 3','n_0 = 4');
title('Repeated Quantization - \gamma_{min} = 10[dB]')


figure;hold all;
n_0 = 1:1:4;
ENR_min = 10; ENR_min_lin = 10^(ENR_min/10);
ENR = 6:0.25:40; ENR_lin = 10.^(ENR/10);
c_opt = ((72/pi^2) * sqrt(2*pi*ENR_min_lin)/ENR_min_lin)^(1/3);
beta = c_opt * exp(ENR_min_lin / 6);

plot(ENR,10*log10(1./(1 + ENR_lin.^2)),'-','LineWidth',4); % desired profile

DistortionBound = zeros(size(ENR));
currBound = zeros(length(n_0),length(ENR));

currIdx = find((0 < ENR_lin) & (ENR_lin < pi^2*beta^2 / 6));
currENR = ENR_lin(currIdx);
DistortionBound(currIdx) = (6./(pi^2 * beta^2 * currENR)) + (1/6)*exp(-currENR/2).*(beta-1)./sqrt(2*pi*currENR);

currIdx = find(((pi^2*beta^2 / 6) < ENR_lin) & (ENR_lin < (pi^2*beta^2 / 6)*(beta^2/18)));
probErr = (1/12)*exp(-ENR_lin(currIdx)/2).*(beta-1)./sqrt(2*pi*ENR_lin(currIdx));
currENR = ENR_lin(currIdx)/18;
DistortionBound(currIdx) = probErr + beta^(-2) * (6./(pi^2 * beta^2 * currENR) + (1/6)*exp(-currENR/2).*(beta-1)./sqrt(2*pi*currENR));

currIdx = find((ENR_lin > (pi^2*beta^2 / 6)*(beta^2/18)));
probErr = (1/12)*exp(-ENR_lin(currIdx)/(2*18)).*(beta-1)./sqrt(2*pi*ENR_lin(currIdx)/18);
currENR = ENR_lin(currIdx)/125;
DistortionBound(currIdx) = probErr + beta^(-4) * (6./(pi^2 * beta^2 * currENR) + (1/6)*exp(-currENR/2).*(beta-1)./sqrt(2*pi*currENR));
% end
plot(ENR,10*log10(DistortionBound),'-*','LineWidth',2); % desired profile

grid on; grid minor;
xlabel('ENR [dB]');
ylabel('Distortion [dB]');
legend('D = 1./(1 + \gamma ^2)','Optimized Scalar Quantization + analog PPM');
title('Repeated Quantization - \gamma_{min} = 10[dB]')