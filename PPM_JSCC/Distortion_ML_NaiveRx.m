close all; clear all; clc;

n = 2:1:10;
P_N0W = 12:0.5:15;
P_N0W_lin = 10.^(P_N0W/10);
tau = 2*pi*exp(1)/12;
currLegend = {};
figure;
for i=1:length(P_N0W)
    currErr = 2.^(-2*n) + (1 - 2.^(-2*n)).*(1 - (1 - qfunc(sqrt(P_N0W_lin(i)/2))).^(n-1));
    currOPTA = tau * (1/(1 + P_N0W_lin(i))).^n;
    
    plot(n,10*log10(1./currErr),'--','LineWidth',2); hold on;
    %     semilogy(n,currOPTA,'--','LineWidth',2); hold on;
    
    currLegend = [currLegend;strcat('P/(N_{0}W) = ',num2str(P_N0W(i)))];
end
grid on; grid minor;
xlabel('n'); ylabel('SDR [dB]');
legend(currLegend);
title('JSCC With PPM, uniform source');

% semilogy(n,currOPTA,'--','LineWidth',2); hold on;

n = 2:1:8;
P_N0W = 5:0.25:20;
P_N0W_lin = 10.^(P_N0W/10);
currLegend = {};
figure;
for i=1:length(n)
    currErr = 2.^(-2*n(i)) + (1 - 2.^(-2*n(i))).*(1 - (1 - qfunc(sqrt(P_N0W_lin/2))).^(n(i)-1));
    
    plot(P_N0W,10*log10(currErr),'--','LineWidth',2); hold on;
    %     semilogy(n,currOPTA,'--','LineWidth',2); hold on;
    
    currLegend = [currLegend;strcat('n = ',num2str(n(i)))];
end
grid on; grid minor;
xlabel('P/(N_{0}W)'); ylabel('Distortion [dB]');
legend(currLegend);
title('JSCC With PPM, uniform source');
hold on; plot(P_N0W,10*log10(tau * exp(-1*P_N0W_lin)),'-*');
ylim([-50 -5])

