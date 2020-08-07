
clear all; clc; close all

N = 1;
P = 1;
rho = 0.9;
snr = 7:1:20;


delta = 2.225; % 2.25
beta = 1.07; % 1.05

delta_t = 2.25;
beta_t = 1.05;

Niters = 2*1e4;

SDR_Kochman = zeros(size(snr));
SDR_Tent = zeros(size(snr));
SDR_Linear = zeros(size(snr));
modErrFinal = zeros(size(snr));

snrLin = 10.^(snr/10);
for i=1:length(snr)
    
%     rho = snrLin(i)/(1 + snrLin(i));
    
    curr_SDR_Linear = 0;
    curr_SDR_Kochman = 0;
    curr_SDR_Tent = 0;

    kochman_err = zeros(1,Niters);kochmanModuloErr = zeros(1,Niters);
    tent_err = zeros(1,Niters);
    for n=1:Niters
        
        %% Generate source + side information
        x = sqrt(P)*randn(1,N);
        n_tilde = sqrt(P - rho^2)*randn(1,N);
        y = rho*x + n_tilde;
        w = (1/sqrt(snrLin(i)))*randn(size(x));
        sigmaW = sqrt(1/snrLin(i));
        
        
        %% Encoding + Channel
        % linear
        RxSig_linear = x + w;
        
        % Kochman Zamir (simple modulo)
        TxSig = mod(beta*x + delta,2*delta) - delta;
        RxSig_kochman = TxSig + w;
        
        % Tent Mapping
        TxSig = mod(beta_t*x + delta_t,2*delta_t) - delta_t;
        tent_idx = find(mod(round((beta_t*x - TxSig)/(2*delta_t)),2) == 1);
        TxSig(tent_idx) = -1*TxSig(tent_idx);
        RxSig_tent = TxSig + w;

        %% Decoding
        %% linear
        x_hat = ((rho/snrLin(i))*y + (1-rho^2)*RxSig_linear)/(1 + 1/snrLin(i) - rho^2);
        curr_SDR_Linear = curr_SDR_Linear + sum((x_hat - x).^2);
        
        %% Kochman
        [x_hat_kochman,modError] = Kochman_Decoder(RxSig_kochman,y,beta,delta,rho,sigmaW);
        curr_SDR_Kochman = curr_SDR_Kochman + sum((x_hat_kochman - x).^2);
        kochman_err(n) = x_hat_kochman - x;
        kochmanModuloErr(n) = modError;
        
        %% Tent 
        x_hat_tent = FullMMSE_decoder_tent(RxSig_tent,y,beta,delta,rho,sigmaW);
        curr_SDR_Tent = curr_SDR_Tent + sum((x_hat_tent - x).^2);
        tent_err(n) = x_hat_tent - x;
    end
    
    SDR_Linear_Final = curr_SDR_Linear/(N*Niters);
    SDR_Linear(i) = 10*log10(1/SDR_Linear_Final);
    
    SDR_Tent_Final = curr_SDR_Tent/(N*Niters);
    SDR_Tent(i) = 10*log10(1/SDR_Tent_Final);

    SDR_Kochman_Final = curr_SDR_Kochman/(N*Niters);
    SDR_Kochman(i) = 10*log10(1/SDR_Kochman_Final);
    
    modErrFinal(i) = sum(kochmanModuloErr)/(N*Niters);
    
    display(strcat('Finished SNR = ',num2str(snr(i))));

end
figure;hold all
plot(snr , 10*log10((1 + snrLin)/(1 - rho^2)),'LineWidth',2)
% plot(snr , SDR_Tent,'b-*','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
plot(snr , SDR_Kochman,'--gs','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
plot(snr , 10*log10((1 + (1-rho^2)*snrLin)/(1-rho^2)),'r--','LineWidth',2)

grid on; grid minor;
xlabel('SNR [dB]'); ylabel('SDR [dB]');
legend('SDR OPTA - SI case : (1 + SNR)/(1-\rho^2)',...
    'Kochman',...
    'Linear Scheme');
title(strcat('Distortion Vs SNR, \rho = ',num2str(rho)));
ylim([10 20]);

figure;semilogy(snr,modErrFinal,'LineWidth',2);
grid on; grid minor;
xlabel('snr [dB]'); ylabel('Prob');
title('Modulo Error Probability');


function [s_hat,modError] = Kochman_Decoder(RxSig,y,beta,delta,rho,sigmaW)

modError = 0;
if abs(RxSig - beta*rho*y) > delta
    modError = 1;
end
mod_sig = mod(RxSig - beta*rho*y + delta,2*delta) - delta;
alpha_s = (beta*(1-rho^2)) / (beta^2 * (1-rho^2) + sigmaW^2);
s_hat = rho*y + alpha_s * mod_sig;

end

function [X_hat] = FullMMSE_decoder_tent(RxSig,y,beta,delta,rho,sigmaW)

dx = 0.001;
x = (-8:dx:8);

% start with 2D calculations
f_x = mod(beta*x + delta,2*delta) - delta;
tent_idx = find(mod(round((beta*x - f_x)/(2*delta)),2) == 1);
f_x(tent_idx) = -1*f_x(tent_idx);


expTerm1 = -1*(0.5*(x).^2);
expTerm2 = - (0.5/(1-rho^2))*(y - rho*x).^2;
expTerm3 = - 0.5*(1/sigmaW^2)*(RxSig - f_x).^2;
Integrand = exp(expTerm1 + expTerm2 + expTerm3);

num = sum(Integrand.*x,2);
den = sum(Integrand);
X_hat = num/den;

end

