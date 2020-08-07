% close all;
clear all; clc; close all

%% lambda = arbitrary < 1
sigma_v = 1;
lambda = 0.9;
sigma_w = 1/4;

d = 0.0001:0.000025:sigma_w^2;
d_vic = 0.0000125:0.000025:sigma_w^2;
% R Victoria
R_vic = 0.5 * log(1./(d_vic.*(1/sigma_w^2 + 1./(lambda^2*d_vic + sigma_v^2))));
R_vic(R_vic <= 0) = 0;

R_noSI = 0.5 * log(lambda^2 + sigma_v^2./d);
R_noSI(R_noSI <= 0) = 0;

d_new = 0.0001:0.00001:(sigma_v^2 * sigma_w^2)/(sigma_v^2 + sigma_w^2);
tempR = zeros(length(d_new),1);
for i=1:length(d_new)
    
    if i == length(d_new)
        disp('Fuck My Life');
    end
    sigma_z_sigma_n = (sigma_v^2 * d_new(i) + lambda^2 * (d_new(i))^2)/(sigma_v^2 - (1 - lambda^2)*d_new(i));
    sigma_z_2 = sigma_w^2 * sigma_z_sigma_n / (sigma_w^2 - sigma_z_sigma_n);
    
    p = [lambda^2,(sigma_v^2 + (1 - lambda^2)*sigma_z_2),-1*sigma_v^2 * sigma_z_2];
    r = roots(p);
    tempR(i) = max(0,0.5 * log(lambda^2 + sigma_v^2 / r(r>0)));
end
R_elGamal = tempR;
R_elGamal(end) = 0;

% convex hull of R(D|U)
currR = [50 reshape(R_elGamal,1,[])];
currD = [50 reshape(d_new,1,[])];
convHull_Idx = convhull(currD,currR);


figure;hold all
plot(d,R_noSI,'-.','LineWidth',2)
plot(d_new,R_elGamal,':','LineWidth',2)
plot(currD(convHull_Idx(:).'),currR(convHull_Idx(:).'),'--','LineWidth',2.5)
plot(d_vic,R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('D [Linear]'); ylabel('Rate [bits]');
title(strcat('Causal Rate Distortion, \lambda = ',num2str(lambda)));
legend('No SI','Rx SI: Gaussian Test Channel','Rx SI: Convexified','Two-Sided SI');
ylim([0 3.35]); xlim([0 max(d_new)])

figure;hold all
plot(10*log10(d),R_noSI,'-.','LineWidth',2)
plot(10*log10(d_new),R_elGamal,':','LineWidth',2)
plot(10*log10(currD(convHull_Idx(:).')),currR(convHull_Idx(:).'),'--','LineWidth',2.5)
plot(10*log10(d_vic),R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('D [dB]'); ylabel('Rate [bits]');
title(strcat('Causal Rate Distortion, \lambda = ',num2str(lambda)));
legend('No SI','Rx SI: Gaussian Test Channel','Rx SI: Convexified','Two-Sided SI');
ylim([0 3.35]); xlim([-30 10*log10(max(d_new))])


figure;
subplot(121); hold all
plot(d,R_noSI,'-.','LineWidth',2)
plot(d_new,R_elGamal,':','LineWidth',2)
plot(currD(convHull_Idx(:).'),currR(convHull_Idx(:).'),'--','LineWidth',2.5)
plot(d_vic,R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('Average Distortion [Linear]'); ylabel('Average Rate [bits]');
% title(strcat('Causal Rate Distortion, \lambda = ',num2str(lambda)));
legend('No SI','Rx SI: Gaussian Test Channel','Rx SI: Convexified','Two-Sided SI');
ylim([0 3.35]); xlim([0 max(d_new)])



subplot(122);
hold all;
plot(10*log10(d),R_noSI,'-.','LineWidth',2)
plot(10*log10(d_new),R_elGamal,':','LineWidth',2)
plot(10*log10(currD(convHull_Idx(:).')),currR(convHull_Idx(:).'),'--','LineWidth',2.5)
plot(10*log10(d_vic),R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('Average Distortion [dB]'); ylabel('Average Rate [bits]');
% title(strcat('Causal Rate Distortion, \lambda = ',num2str(lambda)));
ylim([0 3.35]); xlim([-40 10*log10(max(d_new))])

    