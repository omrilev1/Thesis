clear all; close all; clc;

snr = 10:1:80;
snrLin = 10.^(snr/10);
rho = 0.9;

s = 0.000125:0.000125:4;

lambda = linspace(0.6,18,1000);
MSEbound = zeros(length(lambda),length(snr));
MSEbound_SI = zeros(length(lambda),length(snr));
MSEbound_SI_nonUniform = zeros(length(lambda),length(snr));

MSE_nominal = zeros(length(lambda),length(snr));
MSE_nominal_SI = zeros(length(lambda),length(snr));
MSE_nominal_SI_nonUniform = zeros(length(lambda),length(snr));

Mc_overM = 0.99; % continuous almost everywhere

m = 2;

% Spiral parameters 
sigma_x = 0.25; % for the spiral calculations, later we'll compenste over this
eta = 0.16; 
Delta = 1e-3:1e-3:1.75; 
MSE_spiral = zeros(1,length(snr));
bestDelta = zeros(1,length(snr));

for i=1:length(snr)
    
    anomProb = exp(-2*lambda);
    r_th = (1/sqrt(snrLin(i))) * qfuncinv(0.5*anomProb);

%     R = log(pi*(r_th + 1).^2 ./ r_th);
    R = log(pi*(r_th + 1).^2 ./ r_th);

%     theta = 0:0.01:2*pi; theta = reshape(theta,[],1);
%     covering_radius = pi*(1 + abs(cos(theta))*r_th).*(abs(sin(theta))*r_th + rho/sqrt(snrLin(i)*(1-rho^2)));
%     minRadius = min(covering_radius,[],1);
%     R_SI = log(minRadius./r_th);
    R_SI = log(pi*(r_th/sqrt(2) + 1).*(r_th/sqrt(2) + sqrt(rho^2/((1 - rho^2)*snrLin(i)))) ./ r_th);

    
    
    
    firstTerm = repmat(2^(1+m)*(s(:).^m),1,length(lambda)).*(qfunc(sqrt(snrLin(i)) * s(:)) - Mc_overM * anomProb);
    
    currMSE_SI = firstTerm.*exp(-m*repmat(R_SI,length(s),1));
    currMSE = firstTerm.*exp(-m*repmat(R,length(s),1));

    MSEbound_SI(:,i) = max(currMSE_SI,[],1);
    MSEbound(:,i) = max(currMSE,[],1);
    
    MSE_nominal(:,i) = MSEbound(:,i) .* (1 - anomProb(:)) + 1 * anomProb(:);
    MSE_nominal_SI(:,i) = MSEbound_SI(:,i) .* (1 - anomProb(:)) + 1 * anomProb(:);
    
    % Find best Spiral SDR for that case 
    alpha = eta*sqrt(2*pi^5) ./ (Delta*sigma_x*(1 - exp(-1/(2*sigma_x^2))));
    P_r = (1 - erf(Delta/(2*sqrt(2)/sqrt(snrLin(i)))));
    eps_th = (1./(sqrt(pi)*(alpha.^2))) .* P_r .* (erf(sqrt(2)/(2*sigma_x))*(4*sqrt(pi)*sigma_x^2*(alpha.^2) + eta^2*(Delta.^2)*(pi^(4.5))) - ...
        4*sqrt(2)*sigma_x*alpha*exp(-1/(2*sigma_x^2)).*(alpha + 2*eta*Delta*(pi^2)) + 8*sqrt(2)*eta*(pi^2)*sigma_x*Delta.*alpha);
    [MSE_spiral(i),minIdx] = min(eps_th + (1/snrLin(i))./(alpha.^2));
    bestDelta(i) = Delta(minIdx); 

end

% gaussian channel with average power constratint
bestSDR = max(10*log10(1./MSE_nominal),[],1);
bestSDR_SI = max(10*log10(1./MSE_nominal_SI),[],1);
SDR_Spiral = sigma_x^2 ./ MSE_spiral;

% Plot the SI case 
figure;
subplot(121); hold all 
plot(snr,10*log10((1 + snrLin)/(1 - rho^2)),'-*','LineWidth',2);
plot(snr,bestSDR_SI,'-o','LineWidth',2);
plot(snr,min(10*log10((1 + snrLin)/(1 - rho^2)),bestSDR_SI),'-.g','LineWidth',2)
plot(snr,10*log10((1 + (1-rho^2)*snrLin)/(1-rho^2)),'-*','LineWidth',2);
xlabel('snr [dB]'); ylabel('SDR [dB]');
legend('Wyner Ziv n = \infty','New Bound','overall Bound','Linear Scheme');
title('SDR bound for JSCC with Rx SI');
grid on; grid minor;

subplot(122); hold all 
plot(snr(25:end),10*log10((1 + snrLin(25:end))/(1 - rho^2)),'-*','LineWidth',2);
plot(snr(25:end),bestSDR_SI(25:end),'-o','LineWidth',2);
plot(snr(25:end),min(10*log10((1 + snrLin(25:end))/(1 - rho^2)),bestSDR_SI(25:end)),'-.g','LineWidth',2)
plot(snr(25:end),10*log10((1 + (1-rho^2)*snrLin(25:end))/(1-rho^2)),'-*','LineWidth',2);
xlabel('snr [dB]'); ylabel('SDR [dB]');
legend('Wyner Ziv n = \infty','New Bound','overall Bound','Linear Scheme');
title('Zoom-in high SNR region');
grid on; grid minor;

% Plot the general 1:n case 
figure;
subplot(121); hold all 
plot(snr,10*log10((1 + snrLin).^2),'-*','LineWidth',2);
plot(snr,bestSDR,'-o','LineWidth',2);
plot(snr,10*log10(1 + 2*snrLin),'-*','LineWidth',2);
plot(snr,10*log10(SDR_Spiral),'--','LineWidth',2);
xlabel('snr [dB]'); ylabel('SDR [dB]');
legend('SDR OPTA','New Bound','Linear Transmission','Best Archimedean spiral : analytical');
title('SDR bound for 1:2 JSCC');

grid on; grid minor;

subplot(122); hold all 
plot(snr(25:end),10*log10((1 + snrLin(25:end)).^2),'-*','LineWidth',2);
plot(snr(25:end),bestSDR(25:end),'-o','LineWidth',2);
plot(snr(25:end),10*log10(1 + 2*snrLin(25:end)),'-*','LineWidth',2);
plot(snr(25:end),10*log10(SDR_Spiral(25:end)),'--','LineWidth',2);
xlabel('snr [dB]'); ylabel('SDR [dB]');
legend('SDR OPTA','New Bound','Linear Transmission','Best Archimedean spiral : analytical');

title('Zoom-in high SNR region');
grid on; grid minor;
