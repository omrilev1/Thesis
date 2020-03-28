% close all;
clear all; clc;
t = 1000;

%% lambda = 0
sigma_v = 1;
sigma_w = 1/4;

d = 0.0001:0.000025:sigma_w^2;
d_vic = 0.0000125:0.000025:sigma_w^2;
% R Victoria
R_vic = 0.5 * log(1./(d_vic.*(1/sigma_w^2 + 1./sigma_v^2)));
R_vic(R_vic <= 0) = 0;

R_noSI = 0.5 * log(sigma_v^2./d);
R_noSI(R_noSI <= 0) = 0;

% Our R :
dVec = d;
tempR = zeros(length(d),t);
tempR1 = zeros(length(d),t);
tempR2 = zeros(length(d),t);

% Solve for D_c
[minVal,idx_Dc] = min(abs(-(1./(2*d.^2)).*(1./(1./d - log(2)/sigma_w^2)) - 0.5*log(sigma_v^2* (1./d - 1/sigma_w^2)) ./ (d - sigma_v^2*sigma_w^2/(sigma_v^2 + sigma_w^2))));
Dc = d(idx_Dc);

for i=1:(t-1)
    
    % Find D_c
    
    % Regime 1
    tempR(dVec > sigma_v^2 * sigma_w^2 / (sigma_v^2 + sigma_w^2),i) = 0;
    
    % Regime 2 :
    idx2 = find((dVec >=0) & (dVec < Dc));
    tempR(idx2,i) = 0.5*log(sigma_v^2* (1./dVec(idx2) - 1/sigma_w^2));
    
    % Regime 3 :
    idx3 = find((dVec >= Dc) & (dVec < sigma_v^2 * sigma_w^2 / (sigma_v^2 + sigma_w^2)));
    
    tempR(idx3,i) = (dVec(idx3) - sigma_v^2*sigma_w^2/(sigma_v^2 + sigma_w^2)).*(-(1./(2.* (Dc.^2))).*(1./(1./Dc - log(2)/sigma_w^2)));
    
    
end
R_elGamal = sum(tempR,2)/t;


% R bound 
sigma_v_sigma_w = sigma_v^2 * sigma_w^2 / (sigma_v^2 + sigma_w^2);

figure;hold all
plot(d,R_noSI,'-.','LineWidth',2)
plot(d,max(0.5*log(sigma_v^2* (1./d - 1/sigma_w^2)),0),':','LineWidth',2)
plot(d,R_elGamal,'--','LineWidth',2)
plot(d_vic,R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('Average Distortion [Linear]'); ylabel('Average Rate [bits]');
% title('Causal Rate Distortion, \lambda = 0');
legend('Without SI','SI Causally @Rx','SI Causally @Rx - Convexification','Two-Sided SI (Wyner-Ziv)');


figure;hold all
plot(10*log10(d),R_noSI,'-.','LineWidth',2)
plot(10*log10(d),max(0.5*log(sigma_v^2* (1./d - 1/sigma_w^2)),0),':','LineWidth',2)
plot(10*log10(d),R_elGamal,'--','LineWidth',2)
plot(10*log10(d_vic),R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('Average Distortion [dB]'); ylabel('Average Rate [bits]');
% title('Causal Rate Distortion, \lambda = 0');
legend('Without SI','SI Causally @Rx','SI Causally @Rx - Convexification','Two-Sided SI (Wyner-Ziv)');

figure;
subplot(121);
hold all
plot(d,R_noSI,'-.','LineWidth',2)
plot(d,max(0.5*log(sigma_v^2* (1./d - 1/sigma_w^2)),0),':','LineWidth',2)
plot(d,R_elGamal,'--','LineWidth',2)
plot(d_vic,R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('Average Distortion [Linear]'); ylabel('Average Rate [bits]');
% title('Causal Rate Distortion, \lambda = 0');
legend('Without SI','SI Causally @Rx','SI Causally @Rx - Convexification','Two-Sided SI (Wyner-Ziv)');


subplot(122);
hold all;
plot(10*log10(d),R_noSI,'-.','LineWidth',2)
plot(10*log10(d),max(0.5*log(sigma_v^2* (1./d - 1/sigma_w^2)),0),':','LineWidth',2)
plot(10*log10(d),R_elGamal,'--','LineWidth',2)
plot(10*log10(d_vic),R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('Average Distortion [dB]'); ylabel('Average Rate [bits]');
% title('Causal Rate Distortion, \lambda = 0');

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

tempR = zeros(length(d),t);
for i=1:(t-1)
    if i==1
        tempR(:,1) = R_elGamal;
    else
        tempR(:,i) = 0.5*log(sigma_v^2 * (1./d - sigma_v^2./(lambda^2 * d + sigma_v^2) - sigma_v^2/sigma_w^2) + lambda^2*(1 - 2.^(-2*reshape(tempR(:,i-1),1,[]))) + 1);
    end
    tempR(tempR(:,i) <= 0,i) = 0;
end
R_elGamal = sum(tempR,2)/t;

% convex hull of R(D|U)
currR = [50 reshape(R_elGamal,1,[])];
currD = [50 reshape(d,1,[])];
convHull_Idx = convhull(currD,currR);

figure;hold all
plot(d,R_noSI,'-.','LineWidth',2)
plot(d,R_elGamal,':','LineWidth',2)
plot(currD(convHull_Idx(:).'),currR(convHull_Idx(:).'),'--','LineWidth',2.5)
plot(d_vic,R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('D [Linear]'); ylabel('Rate [bits]');
title(strcat('Causal Rate Distortion, \lambda = ',num2str(lambda)));
legend('Without SI','SI Causally @Rx','SI Causally @Rx - Convexification','Two-Sided SI (Wyner-Ziv)');
ylim([0 3.35]); xlim([0 max(d)])

figure;hold all
plot(10*log10(d),R_noSI,'-.','LineWidth',2)
plot(10*log10(d),R_elGamal,':','LineWidth',2)
plot(10*log10(currD(convHull_Idx(:).')),currR(convHull_Idx(:).'),'--','LineWidth',2.5)
plot(10*log10(d_vic),R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('D [dB]'); ylabel('Rate [bits]');
title(strcat('Causal Rate Distortion, \lambda = ',num2str(lambda)));
legend('Without SI','SI Causally @Rx','SI Causally @Rx - Convexification','Two-Sided SI (Wyner-Ziv)');
ylim([0 3.35]); xlim([-30 10*log10(max(d))])


figure;
subplot(121); hold all
plot(d,R_noSI,'-.','LineWidth',2)
plot(d,R_elGamal,':','LineWidth',2)
plot(currD(convHull_Idx(:).'),currR(convHull_Idx(:).'),'--','LineWidth',2.5)
plot(d_vic,R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('Average Distortion [Linear]'); ylabel('Average Rate [bits]');
% title(strcat('Causal Rate Distortion, \lambda = ',num2str(lambda)));
legend('Without SI','SI Causally @Rx','SI Causally @Rx - Convexification','Two-Sided SI (Wyner-Ziv)');
ylim([0 3.35]); xlim([0 max(d)])



subplot(122);
hold all;
plot(10*log10(d),R_noSI,'-.','LineWidth',2)
plot(10*log10(d),R_elGamal,':','LineWidth',2)
plot(10*log10(currD(convHull_Idx(:).')),currR(convHull_Idx(:).'),'--','LineWidth',2.5)
plot(10*log10(d_vic),R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('Average Distortion [dB]'); ylabel('Average Rate [bits]');
% title(strcat('Causal Rate Distortion, \lambda = ',num2str(lambda)));
ylim([0 3.35]); xlim([-40 10*log10(max(d))])
