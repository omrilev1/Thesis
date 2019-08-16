close all; clear all; clc;
% channel with feedback try


snr = [6]; % snr to check
P = 1; % input power constraint
N_feedback = 4:1:24;
delta_snr = 14;

delta_theory = [3.25 5.85 8.7];

dist_modulo = zeros(length(snr),length(N_feedback));
dist_linear = zeros(length(snr),length(N_feedback));
dist_linear_optimized = zeros(length(snr),length(N_feedback));
theory_modulo = zeros(length(snr),length(N_feedback));
theory_linear = zeros(length(snr),length(N_feedback));

N = 1e4;
snrLin = 10.^(snr/10);
snrB_Lin = 10.^((snr + delta_snr)/10);
P_alias = 1e-6;
for i=1:length(snr)
    
    for j=1:length(N_feedback)
        curr_dist_modulo = 0; curr_dist_linear = 0; curr_dist_linear_optimized = 0;
        for n = 1:N
            theta = randn; % gaussian source with unit power constraint
            
            [x_final,MSE_final] = feedback_over_awgn(theta,P,P,snr(i),delta_snr,N_feedback(j),P_alias,'modulo');
            [x_final_linear,MSE_final_linear] = feedback_over_awgn(theta,P,P,snr(i),delta_snr,N_feedback(j),P_alias,'linear');
            [x_final_linear_optimized,MSE_final_linear_optimized] = feedback_over_awgn(theta,P,P,snr(i),delta_snr,N_feedback(j),P_alias,'linear optimized');
            
            curr_dist_modulo           = curr_dist_modulo + abs(x_final-theta).^2 ;
            curr_dist_linear           = curr_dist_linear + abs(x_final_linear-theta).^2;
            curr_dist_linear_optimized = curr_dist_linear_optimized + abs(x_final_linear_optimized-theta).^2;
            
        end
        dist_modulo(i,j) = 10*log10(curr_dist_modulo/N);
        dist_linear(i,j) = 10*log10(curr_dist_linear/N);
        dist_linear_optimized(i,j) = 10*log10(curr_dist_linear_optimized/N);
        theory_modulo(i,j) = 10*log10(MSE_final);
        theory_linear(i,j) = 10*log10(MSE_final_linear);
        
        display(strcat('Finished feedback = ',num2str(N_feedback(j))))
    end
    
    figure; hold all
    plot(N_feedback,dist_modulo(i,:),'-*','LineWidth',2.5)
    plot(N_feedback,dist_linear(i,:),'-o','LineWidth',2.5)
    plot(N_feedback,dist_linear_optimized(i,:),'-gs','LineWidth',2.5)
    plot(N_feedback,10*log10(1./(N_feedback*(snrLin + snrB_Lin))),'-','LineWidth',3.5)
    %     plot(N_feedback,theory_modulo(i,:) - delta_theory(i),'--','LineWidth',2.5)
    %     plot(N_feedback,theory_linear(i,:),'-*','LineWidth',2.5)
    
    curString = strcat('MSE For snr = ',num2str(snr(i)),' [dB]');
    grid on; grid minor;
    title(strcat('SNR = ',num2str(snr(i)),' [dB] \DeltaSNR = ',num2str(delta_snr),'[dB]'))
    legend({'Modulo','Linear feedback of $\hat{\theta}_{n}$','Linear transmission with optimized coefficents',...
        'Linear Transmission theory $\frac{1}{N*(snr_{F} + snr_{B})}$'},'FontSize',16)
    set(legend,'Interpreter','latex')
    
    xlabel('N feedback iterations')
    ylabel('MSE')
    
end

%{
% Feedback of theta estimate
epsilon = zeros(1,20);
epsilon(1) = 1/snrLin;

for i=2:20
    
    epsilon(i) = (epsilon(i-1) + (1 + epsilon(i-1))/snrB_Lin)/snrLin + (1 + epsilon(i-1))/snrB_Lin;
end
figure;hold all; plot(1:20,10*log10(epsilon),'LineWidth',2)
plot(1:20,10*log10(((1/snrLin + 1)/snrB_Lin)/(1 - 1/snrLin - 1/snrB_Lin - 1/(snrLin*snrB_Lin))) * ones(1,20),'LineWidth',2)
grid on; grid minor;
title({'Theoretical MSE for Linear Feedback(using Unbiased estimates)'; strcat('SNR_F = ',num2str(snr),' SNR_B = ',num2str(snr + delta_snr))});
legend('Theoretical Curve (unbiased)','Steady State Theory')


% Feedback of received signal
epsilon2 = zeros(1,20);zeta = zeros(1,20);epsilon_hat = zeros(1,20);
epsilon2(1) = 1/snrLin;
zeta(1) = (1/snrB_Lin)*(1 + 1/snrLin);
epsilon_hat(1) = (1/snrB_Lin)*(1 + 1/snrLin) + 1/snrLin;
for i=2:20
    
    epsilon2(i) = zeta(i-1) + epsilon_hat(i-1)/snrLin;
    zeta(i) = zeta(i-1) + (1/snrB_Lin)*(1 + 1/snrLin) * epsilon_hat(i-1);
    epsilon_hat(i) = epsilon_hat(i-1) * ((1/snrB_Lin)*(1 + 1/snrLin) + 1/snrLin);
end
figure;hold all;
plot(1:20,10*log10(epsilon2),'LineWidth',2)
plot(1:20,10*log10(epsilon),'--','LineWidth',2)
plot(1:20,10*log10(((1/snrLin + 1)/snrB_Lin)/(1 - 1/snrLin - 1/snrB_Lin - 1/(snrLin*snrB_Lin))) * ones(1,20),'LineWidth',2)

grid on; grid minor;
title({'Theoretical MSE for Linear Feedback(using Unbiased estimates)'; strcat('SNR_F = ',num2str(snr),' SNR_B = ',num2str(snr + delta_snr))});
legend('Feedback of \theta \hat','Feedback of \epsilon \hat','Steady State')
%}

