% close all;
clear all; clc;
t = 1000;

%% lambda = 0
sigma_v = 1;
sigma_w = 1/4;

d = 0.001:0.00025:sigma_w^2;
d_vic = 0.000125:0.00025:sigma_w^2;
% R Victoria
R_vic = 0.5 * log(1./(d_vic.*(1/sigma_w^2 + 1./sigma_v^2)));
R_vic(R_vic <= 0) = 0;

R_noSI = 0.5 * log(sigma_v^2./d);
R_noSI(R_noSI <= 0) = 0;

% Our R :
dVec = d;
tempR = zeros(length(d),1);
tempR1 = zeros(length(d),1);
tempR2 = zeros(length(d),1);

% Solve for D_c
[minVal,idx_Dc] = min(abs(-(1./(2*d.^2)).*(1./(1./d - log(2)/sigma_w^2)) - 0.5*log(sigma_v^2* (1./d - 1/sigma_w^2)) ./ (d - sigma_v^2*sigma_w^2/(sigma_v^2 + sigma_w^2))));
Dc = d(idx_Dc);


% Find D_c

% Regime 1
tempR(dVec > sigma_v^2 * sigma_w^2 / (sigma_v^2 + sigma_w^2),1) = 0;

% Regime 2 :
idx2 = find((dVec >=0) & (dVec < Dc));
tempR(idx2,1) = 0.5*log(sigma_v^2* (1./dVec(idx2) - 1/sigma_w^2));

% Regime 3 :
idx3 = find((dVec >= Dc) & (dVec < sigma_v^2 * sigma_w^2 / (sigma_v^2 + sigma_w^2)));

tempR(idx3,1) = (dVec(idx3) - sigma_v^2*sigma_w^2/(sigma_v^2 + sigma_w^2)).*(-(1./(2.* (Dc.^2))).*(1./(1./Dc - log(2)/sigma_w^2)));


R_elGamal = tempR;

%% Time Sharing 
sigma_z = 0.001:0.001:10;
alpha = 0.001:0.001:1;
[mesh_alpha,mesh_z] = meshgrid(alpha,sigma_z);
sigma_x_sigma_z = mesh_z.^2 ./ (1 + mesh_z.^2);
sigma_x_sigma_z_sigma_n = (sigma_w^2 * mesh_z.^2 ./ (mesh_z.^2 + sigma_w^2 * (1 + mesh_z.^2)));

R_TS = zeros(size(d));
for i=1:length(d)
    validIdx = find(abs(d(i) - mesh_alpha .* sigma_x_sigma_z_sigma_n - (1 - mesh_alpha).*sigma_x_sigma_z) < 1e-4);
       
   R_TS(i) = min(max(mesh_alpha(validIdx)/2 .* log(1 + 1./(mesh_z(validIdx)).^2),0));
    
end


% R bound
sigma_v_sigma_w = sigma_v^2 * sigma_w^2 / (sigma_v^2 + sigma_w^2);

figure;hold all
plot(d,R_noSI,'-.','LineWidth',2)
plot(d,max(0.5*log(sigma_v^2* (1./d - 1/sigma_w^2)),0),':','LineWidth',2)
plot(d,R_elGamal,'--','LineWidth',2)
plot(d,R_TS,'-','LineWidth',2)
plot(d_vic,R_vic,'-','LineWidth',2)
grid on; grid minor;
xlabel('Average Distortion [Linear]'); ylabel('Average Rate [bits]');
% title('Causal Rate Distortion, \lambda = 0');
legend('Without SI','SI Causally @Rx','SI Causally @Rx - Convexification','Two-Sided SI (Wyner-Ziv)');