close all; clear all; clc;
alpha = 0:0.01:1;
SIR = [0,-3,3]; SIRlin = 10.^(SIR/10);
SNR = -10:2:10; SNRlin = 10.^(SNR/10);

for i=1:length(SIR)
    for j=1:length(SNR)
        obj = (1+alpha.^2 * (1+1/SIRlin(i)) + alpha.^2 / SNRlin(j))./((1-alpha).^2 + alpha.^2 / SNRlin(j));
        figure;hold all;
        plot(alpha,obj,'-ro')
        plot(SNRlin(j)./(1+SNRlin(j)) * ones(size(alpha)) , linspace(0,max(obj) + 1,length(alpha)),'-k','LineWidth',1.5)
        grid on; grid minor;
        xlabel('\alpha'); ylabel('objective');
        legend('objective function','\alpha_{MMSE}')
        title(strcat('objective for SIR =',num2str(SIR(i)),' SNR =',num2str(SNR(j))))
    end
end