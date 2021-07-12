% plot Results - 1D Modulo 
close all;
 
SNR = -7:0.25:35;  % 8              % Direct channel SNR
snrLin = 10.^(SNR/10);
profile2 = 1./(1 + (snrLin).^2);


load('PPM_Profile2.mat');
totalEnergy_PPM = totalEnergy; 
load('Linear_Profile2.mat');
totalEnergy_Linear = totalEnergy; 

figure;subplot(1,2,1); 
plot(SNR,10*log10(profile),'-.','LineWidth',2);hold on;
plot(SNR,10*log10(D_Linear),'LineWidth',2); hold on;
plot(SNR,10*log10(D_PPM),'--k','LineWidth',2); hold on;
xlabel('$$\tilde{E}/N$$ [dB]','FontSize',14,'Interpreter','Latex'); ylabel('Distortion [dB]','FontSize',14,'Interpreter','Latex');
lgd = legend({'Profile','Linear','PPM'},'FontSize',10,'TextColor','Black','Location','Best');
grid on; grid minor;

subplot(1,2,2); 
semilogy(SNR,totalEnergy_Linear,'LineWidth',2); hold on;
semilogy(SNR,totalEnergy_PPM,'-.k','LineWidth',2);hold on;
xlabel('$$\tilde{E}/N$$ [dB]','FontSize',14,'Interpreter','Latex');
ylabel('Accumulated Energy/$$\tilde{E}$$','FontSize',14,'Interpreter','Latex');
ylim([0.75 2.2])
lgd = legend({'Linear','PPM'},'FontSize',10,'TextColor','Black','Location','Best');
grid on; grid minor;
% sgtitle({'Quadratic Profile - Scalar Simulation'},'FontSize',14);
 
load('PPM_Profile3.mat');
load('Linear_Profile3.mat');
subplot(1,2,2); 
semilogy(SNR,profile,'-.','LineWidth',2.5);hold on;
semilogy(SNR,D_Linear,'LineWidth',2.5); hold on;
semilogy(SNR,D_PPM,'--k','LineWidth',2.5); hold on;
xlabel('$$\tilde{E}/N$$ [dB]','FontSize',14,'Interpreter','Latex');
lgd = legend({'Profile','Scalar','PPM'},'FontSize',10,'TextColor','Black','Location','Best');
grid on; grid minor;
title({'Third Order Profile', strcat('\Delta = ',num2str(Delta), ', \alpha = ',num2str(alpha))},'FontSize',14);



%% Infinite Dimension 
SNR = -7:0.25:35;  % 8              % Direct channel SNR
snrLin = 10.^(SNR/10);
profile = 1./(1 + (snrLin).^2);


load('PPMInfDim_Profile2.mat');
totalEnergy_PPM = totalEnergy; 
load('LinearInfDim_Profile2.mat');
totalEnergy_Linear = totalEnergy; 
load('PpmInfDim_Anal_Profile2.mat');
totalEnergy_PPM_Anal = totalEnergy; 


figure;subplot(1,2,1); 
plot(SNR,10*log10(profile),'-.','LineWidth',2);hold on;
plot(SNR,10*log10(D_Linear),'LineWidth',2); hold on;
plot(-7:0.25:30,10*log10(D_PPM),'--k','LineWidth',2); hold on;
plot(SNR,10*log10(D_PPM_Anal),'LineWidth',2); hold on;
xlabel('$$\tilde{E}/N$$ [dB]','FontSize',14,'Interpreter','Latex'); ylabel('Distortion [dB]','FontSize',14,'Interpreter','Latex');
lgd = legend({'Profile','Linear','PPM Empiric','PPM Analytic Bound'},'FontSize',10,'TextColor','Black','Location','Best');
grid on; grid minor;

subplot(1,2,2); 
semilogy(SNR,totalEnergy_Linear,'LineWidth',2); hold on;
semilogy(SNR,totalEnergy_PPM_Anal,'-.','LineWidth',2);hold on;
semilogy(-7:0.25:30,totalEnergy_PPM,'--','LineWidth',2);hold on;
xlabel('$$\tilde{E}/N$$ [dB]','FontSize',14,'Interpreter','Latex');
ylabel('Accomulated Energy/$$\tilde{E}$$','FontSize',14,'Interpreter','Latex');
ylim([0.75 2.2])
lgd = legend({'Linear','PPM Analytic','PPM Empiric'},'FontSize',10,'TextColor','Black','Location','Best');
grid on; grid minor;

% sgtitle({'Quadratic Profile - Infinite Dimension Simulation'},'FontSize',14);
