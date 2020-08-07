clear all; close all; clc
% Modulo folding noise probability : 1D

snr = 8; % -16:2:16;
snrLin = 10.^(snr/10);

P = 10; Delta = sqrt(12*P);
beta_temp = 1:0.25:4;
Prob = zeros(length(snr),length(beta_temp));
N = 0.5*1e7;
for i=1:length(snr)
    P_Z = P/snrLin(i);
    alpha = snrLin(i)/(snrLin(i) + 1);
    curr_snr = 1/(snrLin(i) + 1);
    beta = curr_snr * beta_temp;

    % different lattice zoom-in factors
    for j=1:length(beta)

        currZ = alpha*sqrt(P_Z)*randn(1,N) + (alpha - 1)*Delta*(rand(1,N) - 0.5);
    
        % calculate probability 
        Prob(i,j) = sum(abs(currZ) > Delta/2) / length(currZ);
        
    end
    disp(strcat('Finished SNR = ',num2str(snr(i))));
end

Prob(Prob == 0) = 0.5/N;

figure;imagesc(Prob);
colorbar
xlabel('snr [dB]'); ylabel('\beta');

% 
% figure;semilogy(snr,Prob,'LineWidth',2);
% grid on; grid minor;
% xlabel('snr [dB]'); ylabel('Probability');