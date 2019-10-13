% Testing several transmission schemes of gaussian source where the
% receiver is equipped with side information
% - Regular Scaling and analog transmission
% - Tuncel Scheme, where we send the quantized value + quantization error,
% and use MMSE decoder in the reciever
% - Tuncel scheme, where the receiver knows T in a "genie aided" manner.
% (Upper bound on the optimal performance)

% We also plot the results of linear scheme with Tx Side information -
% acheieves the SDR OPTA  = (1 + SNR/(1-rho^2))

clear all; clc;

% uniform quantizer
delta = 3.15;% 1.5/sqrt(2);
type = 'centroid';
partition = (-6.5:1:6.5)*delta;
switch type
    case 'regular'
        codebook = (-7:1:7)*delta;
    case 'centroid'
        [codebook,~] = generateGaussianCodeBook(partition);
        %probs = [logspace(log(1e-3),log(1.15),(length(codebook) - 1)/2) 1 fliplr(logspace(log(1e-3),log(1.15),(length(codebook) - 1)/2))];
        probs = [logspace(log(1e-3),log(0.67),(length(codebook) - 1)/2) 1.035 fliplr(logspace(log(1e-3),log(0.67),(length(codebook) - 1)/2))];
end

N = 1;
P = 1;
rho = 0.9;
snr = 14;
alpha = 0.8;% 0.8;
beta = -1.15;% -1.15;

% improved tuncel
alpha_tilde = 0.75;% 0.75;
beta_tilde = -1.15;% -1.155;

delta_kochman = 2.25;
beta_kochman = 1.05;

centers = [partition(1) (partition(1:end-1) + partition(2:end))/2 partition(end)];
Niters = 8*1e4;

SDR_Tuncel = zeros(size(snr));
SDR_Kochman = zeros(size(snr));
SDR_Linear = zeros(size(snr));
SDR_TuncelGenie = zeros(size(snr));
SDR_TuncelImproved = zeros(size(snr));
SDR_TuncelImproved_Genie = zeros(size(snr));
SDR_TuncelImproved_FullMMSE = zeros(size(snr));
SDR_Tuncel_FullMMSE = zeros(size(snr));

MSE_T_Tuncel = zeros(size(snr));
MSE_T_TuncelImp = zeros(size(snr));
MSE_S_Tuncel = zeros(size(snr));
MSE_S_TuncelImp = zeros(size(snr));

Bias_Tuncel = zeros(size(snr));
Rho_Tuncel = zeros(size(snr));

snrLin = 10.^(snr/10);
for i=1:length(snr)
    
    curr_SDR_Tuncel = 0;
    curr_SDR_Linear = 0;
    curr_SDR_TuncelGenie = 0;
    curr_SDR_Kochman = 0;
    curr_SDR_Tuncel_Improved = 0;
    curr_SDR_Tuncel_Improved_Genie = 0;
    curr_SDR_Tuncel_Improved_FullMMSE = 0;
    curr_SDR_Tuncel_FullMMSE = 0;
    
    curr_Bias_Tuncel = 0;
    curr_rho_Tuncel = 0;
    curr_Bias_Tuncel_Improved = 0;
    curr_rho_Tuncel_Improved = 0;
    Tdist = zeros(1,Niters); Sdist = zeros(1,Niters);
    
    errTuncel = zeros(1,Niters);
    errGenie = zeros(1,Niters);
    errTuncelImp = zeros(1,Niters);
    errTuncelImpGenie = zeros(1,Niters);
    
    Terr_Tuncel = zeros(1,Niters);
    Terr_TuncelImp = zeros(1,Niters);
    
    Serr_Tuncel = zeros(1,Niters);
    Serr_TuncelImp = zeros(1,Niters);
    
    % solve current optimization problem
%         [probsOptim,cost] = findBestBeta(alpha_tilde,delta,codebook(6:10),rho,10^(snr(i)/10));
%         probs(6:10) = -1*probsOptim;
%         beta_tilde = 1;
    for n=1:Niters
        
        %% Generate source + side information
        x = sqrt(P)*randn(1,N);
        n_tilde = sqrt(P - rho^2)*randn(1,N);
        y = rho*x + n_tilde;
        w = (1/sqrt(snrLin(i)))*randn(size(x));
        %% Encoding + Channel
        % linear
        RxSig_linear = x + w;
        
        % Tuncel
        currDistance = abs(x - centers(:));
        [~,minIdx] = min(currDistance,[],1);
        T = codebook(minIdx);
        S = x - T;
        
        TxSig = alpha*T + beta*S;
        RxSig_tuncel = TxSig + w;
        
        % Debug - Manual Power Calibration
        %        currenPower = mean(TxSig.^2);
        %         figure;
        %         subplot(221); histogram(S); title(strcat('S part histogram, E(s^2) = ',num2str(mean(S.^2))));grid on; grid minor;
        %         subplot(222); histogram(TxSig); title(strcat('TxSig histogram, E(TxSig^2) = ',num2str(mean(TxSig.^2))));grid on; grid minor;
        %         subplot(223); scatter(rho*x(1:10:end),TxSig(1:10:end)); title('[TxSig,rho*x]');xlabel('rho*x');ylabel('TxSig');grid on; grid minor;
        %         subplot(224); scatter(y(1:10:end),RxSig_tuncel(1:10:end)); title(strcat('[RxSig,Y], SNR = ',num2str(snr(i)),'[dB]'));xlabel('y');ylabel('RxSig');grid on; grid minor;
        
        % Tuncel Improved
        TxSig = alpha_tilde*T + beta_tilde*probs(minIdx).*S;
        RxSig_tuncelImproved = TxSig + w;
        
        % Debug - Manual Power Calibration
        currenPowerImp = mean(TxSig.^2);
        %                 figure;
        %                 subplot(221); histogram(S); title(strcat('S part histogram, E(s^2) = ',num2str(mean(S.^2))));grid on; grid minor;
        %                 subplot(222); histogram(TxSig); title(strcat('TxSig histogram, E(TxSig^2) = ',num2str(mean(TxSig.^2))));grid on; grid minor;
        %                 subplot(223); scatter(rho*x(1:100:end),TxSig(1:100:end)); title('[TxSig,rho*x]');xlabel('rho*x');ylabel('TxSig');grid on; grid minor;
        %                 subplot(224); scatter(y(1:100:end),RxSig_tuncel(1:100:end)); title(strcat('[RxSig,Y], SNR = ',num2str(snr(i)),'[dB]'));xlabel('y');ylabel('RxSig');grid on; grid minor;
        %
        
        % Kochman Zamir
        curr_dither = 0; % 2*delta_kochman*rand - delta_kochman;
        TxSig = mod(beta_kochman*x + delta_kochman + curr_dither,2*delta_kochman) - delta_kochman;
        RxSig_kochman = TxSig + w;
        
        % Debug - Manual Power Calibration
        %         currenPower = mean(TxSig.^2);
        %         figure;
        %         subplot(211); histogram(TxSig); title(strcat('TxSig histogram, E(TxSig^2) = ',num2str(mean(TxSig.^2))));grid on; grid minor;
        %         subplot(212); scatter(rho*x(1:10:end),TxSig(1:10:end)); title('[TxSig,rho*x]');xlabel('rho*x');ylabel('TxSig');grid on; grid minor;
        %         subplot(213); scatter(y(1:10:end),RxSig_kochman(1:10:end)); title(strcat('[RxSig,Y], SNR = ',num2str(snr(i)),'[dB]'));xlabel('y');ylabel('RxSig');grid on; grid minor;
        %
        
        %% Decoding
        %% linear
        x_hat = ((rho/snrLin(i))*y + (1-rho^2)*RxSig_linear)/(1 + 1/snrLin(i) - rho^2);
        curr_SDR_Linear = curr_SDR_Linear + sum((x_hat - x).^2);
        
        %% Tuncel
        sigmaW = sqrt(1/snrLin(i));
        T_hat = MMSE_decoder_Int(RxSig_tuncel,y,delta,codebook,alpha,beta,rho,sigmaW,NaN);
        Terr_Tuncel(n) = (T_hat - T);
        
        S_hat = MMSE_decoder_Frac(RxSig_tuncel,y,delta,alpha,beta,rho,sigmaW,codebook,NaN);
        x_hat_Tuncel = T_hat + S_hat;
        
        Serr_Tuncel(n) = (S_hat - S);
        
        curr_SDR_Tuncel = curr_SDR_Tuncel + sum((x_hat_Tuncel - x).^2);
        curr_Bias_Tuncel = curr_Bias_Tuncel + x_hat_Tuncel - x ;
        curr_rho_Tuncel = curr_rho_Tuncel + x*(x_hat_Tuncel - x);

        %% Tuncel Full MMSE
        sigmaW = sqrt(1/snrLin(i));
        x_hat_Tuncel = FullMMSE_decoder(RxSig_tuncel,y,alpha,beta,rho,sigmaW,codebook,centers,ones(size(probs)));

        curr_SDR_Tuncel_FullMMSE = curr_SDR_Tuncel_FullMMSE + sum((x_hat_Tuncel - x).^2);
%         curr_Bias_Tuncel_FullMMSE = curr_Bias_Tuncel_FullMMSE + x_hat_Tuncel - x ;
%         curr_rho_Tuncel_FullMMSE = curr_rho_Tuncel_FullMMSE + x*(x_hat_Tuncel - x);
        
        %% Tuncel Genie - the receiver knows T
        S_hat = MMSE_decoder_Frac(RxSig_tuncel,y,delta,alpha,beta,rho,sigmaW,codebook,T);
        x_hat_Tuncel_Genie = T + S_hat;
        %         Serr_Tuncel(n) = S_hat - S;
        
        curr_SDR_TuncelGenie = curr_SDR_TuncelGenie + sum((x_hat_Tuncel_Genie - x).^2);
        
        %% Tuncel Improved
        sigmaW = sqrt(1/snrLin(i));
        T_hat_Improved = MMSE_decoder_Int_Improved(RxSig_tuncelImproved,y,delta,codebook,probs,alpha_tilde,beta_tilde,rho,sigmaW);
        Terr_TuncelImp(n) = T_hat_Improved - T;
        
        S_hat_Improved = MMSE_decoder_Frac_Improved(RxSig_tuncelImproved,y,delta,alpha_tilde,beta_tilde,rho,sigmaW,codebook,probs,NaN);
        x_hat_Tuncel_Improved = T_hat_Improved + S_hat_Improved;
        Serr_TuncelImp(n) = (S_hat_Improved - S);
        
        curr_SDR_Tuncel_Improved = curr_SDR_Tuncel_Improved + sum((x_hat_Tuncel_Improved - x).^2);
        curr_Bias_Tuncel_Improved = curr_Bias_Tuncel_Improved + x_hat_Tuncel_Improved - x ;
        curr_rho_Tuncel_Improved = curr_rho_Tuncel_Improved + x*(x_hat_Tuncel_Improved - x);
        
        %% Tuncel Improved Genie - the receiver knows T
        S_hat_Improved_Genie = MMSE_decoder_Frac_Improved(RxSig_tuncelImproved,y,delta,alpha_tilde,beta_tilde,rho,sigmaW,codebook,probs,T);
        x_hat_Tuncel_Improved_Genie = T + S_hat_Improved;
        %         Serr_TuncelImp(n) = S_hat_Improved_Genie - S;
        
        curr_SDR_Tuncel_Improved_Genie = curr_SDR_Tuncel_Improved_Genie + sum((x_hat_Tuncel_Improved_Genie - x).^2);
        
        %% Tuncel Improved - Full MMSE
        [X_hat_TuncelImproved_FullMMSE] = FullMMSE_decoder(RxSig_tuncelImproved,y,alpha_tilde,beta_tilde,rho,sigmaW,codebook,centers,probs);
        curr_SDR_Tuncel_Improved_FullMMSE = curr_SDR_Tuncel_Improved_FullMMSE + sum((X_hat_TuncelImproved_FullMMSE - x).^2);
        
        
        %% Kochman
        x_hat_kochman = Kochman_Decoder(RxSig_kochman,y,beta_kochman,delta_kochman,rho,sigmaW,curr_dither);
        
        curr_SDR_Kochman = curr_SDR_Kochman + sum((x_hat_kochman - x).^2);
        
        %% statistics
        errTuncel(n) = x_hat_Tuncel - x;
        errGenie(n) = x_hat_Tuncel_Genie - x;
        errTuncelImp(n) = X_hat_TuncelImproved_FullMMSE - x;
        errTuncelImpGenie(n) = x_hat_Tuncel_Improved_Genie - x;
        
        
        if abs(errTuncel(n)) > 1
            disp(strcat('SI noise : ',num2str(n_tilde/sqrt(P - rho^2))));
            disp(strcat('channel noise : ',num2str(w/(1/sqrt(snrLin(i))))));
        end
    end
    SDR_Tuncel_Final = curr_SDR_Tuncel/(N*Niters);
    SDR_Tuncel(i) = 10*log10(1/SDR_Tuncel_Final);
    Bias_Tuncel(i) = curr_Bias_Tuncel/(N*Niters);
    Rho_Tuncel(i) = curr_rho_Tuncel/(N*Niters);
    
    SDR_Linear_Final = curr_SDR_Linear/(N*Niters);
    SDR_Linear(i) = 10*log10(1/SDR_Linear_Final);
    
    SDR_Kochman_Final = curr_SDR_Kochman/(N*Niters);
    SDR_Kochman(i) = 10*log10(1/SDR_Kochman_Final);
    
    SDR_TuncelGenie_Final = curr_SDR_TuncelGenie/(N*Niters);
    SDR_TuncelGenie(i) = 10*log10(1/SDR_TuncelGenie_Final);
    
    SDR_TuncelImproved_Final = curr_SDR_Tuncel_Improved/(N*Niters);
    SDR_TuncelImproved(i) = 10*log10(1/SDR_TuncelImproved_Final);
    
    SDR_TuncelImproved_Final_FullMMSE = curr_SDR_Tuncel_Improved_FullMMSE/(N*Niters);
    SDR_TuncelImproved_FullMMSE(i) = 10*log10(1/SDR_TuncelImproved_Final_FullMMSE);

    SDR_TuncelImproved_Final_Genie = curr_SDR_Tuncel_Improved_Genie/(N*Niters);
    SDR_TuncelImproved_Genie(i) = 10*log10(1/SDR_TuncelImproved_Final_Genie);
    
    SDR_Tuncel_FullMMSE_Final = curr_SDR_Tuncel_FullMMSE/(N*Niters);
    SDR_Tuncel_FullMMSE(i) = 10*log10(1/SDR_Tuncel_FullMMSE_Final);    
    
    
    display(strcat('Finished SNR = ',num2str(snr(i))));
    
    % Errors in T and S
    MSE_T_Tuncel(i) = mean(Terr_Tuncel.^2);
    MSE_T_TuncelImp(i) = mean(Terr_TuncelImp.^2);
    MSE_S_Tuncel(i) = mean(Serr_Tuncel.^2);
    MSE_S_TuncelImp(i) = mean(Serr_TuncelImp.^2);
    
end
figure;hold all
plot(snr , 10*log10((1 + snrLin)/(1 - rho^2)),'LineWidth',2)
plot(snr , SDR_TuncelGenie,'g.-','LineWidth',2)
plot(snr , SDR_Tuncel_FullMMSE,'-kd','LineWidth',2)
plot(snr , SDR_TuncelImproved,'b-*','LineWidth',2)
plot(snr , SDR_TuncelImproved_Genie,'--gs','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
plot(snr , SDR_TuncelImproved_FullMMSE,'--','LineWidth',2)

% plot(snr , SDR_Kochman,'--gs','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
plot(snr , 10*log10((1 + (1-rho^2)*snrLin)/(1-rho^2)),'r--','LineWidth',2)
% plot(snr , 10*log10((1 + (1-rho^2)*beta^2*snrLin)/(1-rho^2)),'--','LineWidth',2)

grid on; grid minor;
xlabel('SNR [dB]'); ylabel('SDR [dB]');
legend('SDR OPTA - SI case : (1 + SNR)/(1-\rho^2)',...
    strcat('Tuncel (MMSE Decoder) With Genie : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha),' \beta = ',num2str(beta)),...
    strcat('Tuncel (MMSE Decoder) : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha),' \beta = ',num2str(beta)),...
    strcat('Tuncel Improved (MMSE Decoder) : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha_tilde),' \beta = ',num2str(beta_tilde)),...
    strcat('Tuncel Improved with Genie : : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha_tilde),' \beta = ',num2str(beta_tilde)),...
    'Linear Scheme','Tuncel Improved with Full MMSE decoder');
title(strcat('Distortion Vs SNR, \rho = ',num2str(rho)));
ylim([10 20]);

figure;hold all
plot(snr , SDR_TuncelGenie,'g.-','LineWidth',2)
plot(snr , SDR_Tuncel,'-kd','LineWidth',2)
plot(snr , SDR_TuncelImproved,'b-*','LineWidth',2)
plot(snr , SDR_TuncelImproved_Genie,'--gs','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
plot(snr , SDR_TuncelImproved_FullMMSE,'--','LineWidth',2)

grid on; grid minor;
xlabel('SNR [dB]'); ylabel('SDR [dB]');
legend(strcat('Tuncel (MMSE Decoder) With Genie : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha),' \beta = ',num2str(beta)),...
    strcat('Tuncel (MMSE Decoder) : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha),' \beta = ',num2str(beta)),...
    strcat('Tuncel Improved (MMSE Decoder) : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha_tilde),' \beta = ',num2str(beta_tilde)),...
    strcat('Tuncel Improved with Genie : : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha_tilde),' \beta = ',num2str(beta_tilde)),...
    'Tuncel Improved with Full MMSE Decoder');

figure;hold all
plot(snr , SDR_Tuncel_FullMMSE,'b-*','LineWidth',2)
plot(snr , SDR_TuncelImproved_FullMMSE,'--gs','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
grid on; grid minor;
xlabel('SNR [dB]'); ylabel('MSE [dB]');
legend('Tuncel with MMSE decoder','Tuncel improved with Full MMSE decoder');

% % Errors in T and S seperately
% figure;hold all
% plot(snr , 10*log10(MSE_T_Tuncel),'LineWidth',2)
% plot(snr , 10*log10(MSE_T_TuncelImp),'g.-','LineWidth',2)
% plot(snr , 10*log10(MSE_S_Tuncel),'b-*','LineWidth',2)
% plot(snr , 10*log10(MSE_S_TuncelImp),'--gs','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
% grid on; grid minor;
% xlabel('SNR [dB]'); ylabel('MSE [dB]');
% legend('Terr - Tuncel','Terr - Tuncel improved','Serr Tuncel','Serr Tuncel Imp');

% figure;hold all
% plot(snr , -1*10*log10(MSE_T_TuncelImp + MSE_S_TuncelImp),'LineWidth',2)
% plot(snr , -1*10*log10(MSE_T_Tuncel + MSE_S_Tuncel),'LineWidth',2)
% plot(snr , SDR_Tuncel,'b-*','LineWidth',2)
% plot(snr , SDR_TuncelImproved,'--gs','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
% plot(snr , SDR_TuncelImproved_FullMMSE,'--','LineWidth',2)
% grid on; grid minor;
% xlabel('SNR [dB]'); ylabel('MSE [dB]');
% legend('Sum of MSEs Tuncel Imroved','Sum of MSEs Tuncel','Tuncel MSE','Tuncel Imp MSE','Tuncel Improved Full MMSE');

function [T_hat] = MMSE_decoder_Int(RxSig_tuncel,y,delta,codebook,alpha,beta,rho,sigmaW,S)

if isnan(S)
    ds = 0.001;
    s = -delta/2:ds:delta/2;
    
    % start with 2D calculations
    calc1 = codebook(:) + s;
    calc2 = alpha*codebook(:) + beta*s;
    
    expTerm1 = -1*(0.5*(calc1).^2);
    expTerm2 = - (0.5/(1-rho^2))*(y - rho*calc1).^2;
    expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;
    
    Integrand = exp(expTerm1 + expTerm2 - expTerm3);
    
    num = sum(Integrand,2);
    T_hat = sum(codebook(:).*num(:),1)./sum(num(:));
else
    % start with 2D calculations
    calc1 = codebook(:) + S;
    calc2 = alpha*codebook(:) + beta*S;
    
    expTerm1 = -1*(0.5*(calc1).^2);
    expTerm2 = - (0.5/(1-rho^2))*(y - rho*calc1).^2;
    expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;
    
    Integrand = exp(expTerm1 + expTerm2 - expTerm3);
    
    num = Integrand;
    T_hat = sum(codebook(:).*num(:),1)./sum(num(:));
end
end

function [T_hat] = MMSE_decoder_Int_Improved(RxSig_tuncel,y,delta,codebook,probs,alpha,beta,rho,sigmaW)

ds = 0.001;
s = -delta/2:ds:delta/2;

% start with 2D calculations
calc1 = codebook(:) + s;
calc2 = alpha*repmat(codebook(:),1,length(s)) + beta*repmat(probs(:),1,length(s)).*repmat(s,length(codebook),1);

expTerm1 = -1*(0.5*(calc1).^2);
expTerm2 = - (0.5/(1-rho^2))*(y - rho*calc1).^2;
expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;

Integrand = exp(expTerm1 + expTerm2 - expTerm3);

num = sum(Integrand,2);
T_hat = sum(codebook(:).*num(:),1)./sum(num(:));

end

function [S_hat] = MMSE_decoder_Frac(RxSig_tuncel,y,delta,alpha,beta,rho,sigmaW,codebook,T)


if isnan(T)
    ds = 0.001;
    s = ones(length(codebook),1) * (-delta/2:ds:delta/2);
    
    % start with 2D calculations
    calc1 = codebook(:) + s;
    calc2 = alpha*codebook(:) + beta*s;
    
    expTerm1 = -1*(0.5*(calc1).^2);
    expTerm2 = - (0.5/(1-rho^2))*(y - rho*calc1).^2;
    expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;
    Integrand = exp(expTerm1 + expTerm2 - expTerm3);
    
    num = sum(sum(Integrand.*s,2),1);
    den = sum(sum(Integrand));
    S_hat = num/den;
else
    ds = 0.001;
    s = (-delta/2:ds:delta/2);
    
    % start with 2D calculations
    calc1 = T + s;
    calc2 = alpha*T + beta*s;
    
    expTerm1 = -1*(0.5*(calc1).^2);
    expTerm2 = - (0.5/(1-rho^2))*(y - rho*calc1).^2;
    expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;
    Integrand = exp(expTerm1 + expTerm2 - expTerm3);
    
    num = sum(Integrand.*s,2);
    den = sum(Integrand);
    S_hat = num/den;
end
end

function [S_hat] = MMSE_decoder_Frac_Improved(RxSig_tuncel,y,delta,alpha,beta,rho,sigmaW,codebook,probs,T)


if isnan(T)
    ds = 0.001;
    s = ones(length(codebook),1) * (-delta/2:ds:delta/2);
    
    % start with 2D calculations
    calc1 = repmat(codebook(:),1,size(s,2)) + s;
    calc2 = alpha*repmat(codebook(:),1,size(s,2)) + beta*repmat(probs(:),1,length(s)).*s;
    
    expTerm1 = -1*(0.5*(calc1).^2);
    expTerm2 = - (0.5/(1-rho^2))*(y - rho*calc1).^2;
    expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;
    Integrand = exp(expTerm1 + expTerm2 - expTerm3);
    
    num = sum(sum(Integrand.*s,2),1);
    den = sum(sum(Integrand));
    S_hat = num/den;
else
    ds = 0.001;
    s = (-delta/2:ds:delta/2);
    
    %     T_idx = find(codebook == T);
    
    % start with 2D calculations
    calc1 = T + s;
    calc2 = alpha*T + beta*probs(codebook == T).*s;
    
    expTerm1 = -1*(0.5*(calc1).^2);
    expTerm2 = - (0.5/(1-rho^2))*(y - rho*calc1).^2;
    expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;
    Integrand = exp(expTerm1 + expTerm2 - expTerm3);
    
    num = sum(Integrand.*s,2);
    den = sum(Integrand);
    S_hat = num/den;
end
end

function [s_hat] = Kochman_Decoder(RxSig,y,beta,delta,rho,sigmaW,curr_dither)

mod_sig = mod(RxSig - beta*rho*y - curr_dither + delta,2*delta) - delta;

alpha_s = (beta*(1-rho^2)) / (beta^2 * (1-rho^2) + sigmaW^2);
s_hat = rho*y + alpha_s * mod_sig;

end

function [codebook,probs] = generateGaussianCodeBook(partition)

dx = 0.001;
x = (-3 + partition(1)):dx:(partition(end) + 3);

pdf = (1/sqrt(2*pi))*exp(-0.5*x.^2);
codebook = zeros(1,length(partition) - 1);
probs = zeros(1,length(partition) - 1);

for i=1:(length(partition) - 1)
    
    curIdx = find((x >= partition(i)) & (x <= partition(i+1)));
    curPDF = pdf(curIdx) / sum(pdf(curIdx));
    curX = x(curIdx);
    codebook(i) = sum(curX.*curPDF);
    probs(i) = sum(pdf(curIdx))/sum(pdf);
    
end

% generate first and last codewords

% first
curIdx = find(x <= partition(1));
curPDF = pdf(curIdx) / sum(pdf(curIdx));
curX = x(curIdx);
codebook = [sum(curX.*curPDF) codebook];
probs(1) = sum(pdf(curIdx))/sum(pdf);

% last codeword
curIdx = find(x >= partition(end));
curPDF = pdf(curIdx) / sum(pdf(curIdx));
curX = x(curIdx);
codebook = [codebook sum(curX.*curPDF)];
probs(end) = sum(pdf(curIdx))/sum(pdf);

end

function [X_hat] = FullMMSE_decoder(RxSig,y,alpha,beta,rho,sigmaW,codebook,centers,probs)

dx = 0.001;
x = (-6:dx:6);

% start with 2D calculations
currDistance = abs(x - centers(:));
[~,minIdx] = min(currDistance,[],1);
T = codebook(minIdx);
S = x - T;

% Tuncel Improved
f_x = alpha*T + beta*probs(minIdx).*S;

expTerm1 = -1*(0.5*(x).^2);
expTerm2 = - (0.5/(1-rho^2))*(y - rho*x).^2;
expTerm3 = - 0.5*(1/sigmaW^2)*(RxSig - f_x).^2;
Integrand = exp(expTerm1 + expTerm2 + expTerm3);

num = sum(Integrand.*x,2);
den = sum(Integrand);
X_hat = num/den;

end


function [currPower,pdf] = calcQuantizationErrorPower(T,codebook,delta,dx)

if nargin <4
    dx = 0.001;
end
index = find(codebook == T);
index = index - ceil(length(codebook)/2);

x = (index - 0.5)*delta : dx : (index + 0.5)*delta;

pdf = dx*(1/sqrt(2*pi))*exp(-0.5*x.^2);
pdf = pdf/sum(pdf);
currPower = sum(((x - T).^2).*pdf);

end
