clear all; close all; clc;

dt = 1e-4;
alpha = 0:0.01:1;
t = 0:dt:1;

f_alpha = zeros(size(alpha));
g_alpha = zeros(size(alpha));
obj_func = zeros(size(alpha));


snr = -10:1:30;
sir = [0];
snrLin = 10.^(snr/10);
sirLin = 10.^(sir/10);

optAlpha_Tx = zeros(length(sir),length(snr));
optAlpha_Rx = zeros(length(sir),length(snr));

for k=1:length(sir)
    for j=1:length(snr)
        for i=1:length(alpha)
            
            t1 = (alpha(i)/2) : dt : 0.5;
            t2 = (alpha(i)/8) * (1:dt:2).^2; % alpha(i)/2 * (0.5:dt:1).^2;
            f_alpha(i) = 2*sum(H_b(t1)) * dt + alpha(i) * sum(H_b(t2)) * dt;
            
            alpha_Rx =  snrLin(j) * (1 + alpha(i)/sirLin(k)) / (1 + (1 + (1/sirLin(k)))*snrLin(j));
            g_alpha(i) = 0.5*log2(2*pi*exp(1)* ((alpha_Rx)^2 + (alpha_Rx - 1)^2 * snrLin(j) + (alpha_Rx - alpha(i))^2 * snrLin(j)/sirLin(k)));
            
            obj_func(i) = g_alpha(i) + f_alpha(i);
        end
        
        figure;
        plot(alpha,obj_func,'-ks','LineWidth',1.5)
        grid on; grid minor;
        xlabel('\alpha_{Tx}'); ylabel('optimization function');
        title(strcat('Objective function Vs alpha , SNR = ',num2str(snr(j)),'[dB] , SIR = ',num2str(sir(k)),'[dB]'))
        
        [~,opt_idx] = min(obj_func);
        optAlpha_Tx(k,j) = alpha(opt_idx);
        optAlpha_Rx(k,j) = snrLin(j) * (1 + alpha(opt_idx) * (1/sirLin(k))) / (1 + (1 + (1/sirLin(k)))*snrLin(j));
        
    end
    
    figure;hold all
    plot(optAlpha_Tx(k,:),optAlpha_Rx(k,:),'-ks','LineWidth',1.5);
    plot(snrLin./(1+snrLin),snrLin./(1+snrLin),'-ro','LineWidth',1.5);
    grid on; grid minor;
    xlabel('\alpha_{Tx}'); ylabel('\alpha_{Rx}');
    title(strcat('optimal (\alpha_{Tx},\alpha_{Rx}) Pairs , for SIR = ',num2str(sir(k))))
    legend('optimization','\alpha_{Rx} = \alpha_{Tx} = \alpha_{MMSE}')
end
function [Hb] = H_b(x)
Hb = -x.*log2(x) - (1-x).*log2(1-x);
Hb(x == 0) = 0;
end