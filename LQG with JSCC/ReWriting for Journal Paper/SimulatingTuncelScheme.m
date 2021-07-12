% Testing several transmission schemes of gaussian source where the
% receiver is equipped with side information
% - Regular Scaling and analog transmission
% - Tuncel Scheme, where we send the quantized value + quantization error,
% and use MMSE decoder in the reciever
% - Tuncel scheme, where the receiver knows T in a "genie aided" manner.
% (Upper bound on the optimal performance)

% We also plot the results of linear scheme with Tx Side information -
% acheieves the SDR OPTA  = (1 +dbstop SNR/(1-rho^2))

clear all; clc;

% uniform quantizer
delta = 2;% 2;% 3.15;
type = 'centroid';
partition = (-6.5:1:6.5)*delta;
switch type
    case 'regular'
        codebook = (-7:1:7)*delta;
    case 'centroid'
        [codebook,~] = generateGaussianCodeBook(partition);
        %probs = [logspace(log(1e-3),log(1.15),(length(codebook) - 1)/2) 1 fliplr(logspace(log(1e-3),log(1.15),(length(codebook) - 1)/2))];
        probs = [logspace(log(1e-3),log(1.15),(length(codebook) - 1)/2) 1 fliplr(logspace(log(1e-3),log(1.15),(length(codebook) - 1)/2))];
end

N = 1;
P = 1;
rho = 0.95;
snr = -10:1:35;
alpha = 0.85;% 0.8;
beta =  -1.355;% -1.15;

centers = [partition(1) (partition(1:end-1) + partition(2:end))/2 partition(end)];
Niters = 4*1e5;


SDR_Linear = zeros(size(snr));
SDR_Tuncel_FullMMSE = zeros(size(snr));

snrLin = 10.^(snr/10);
for i=1:length(snr)
    
    curr_SDR_Tuncel_FullMMSE = zeros(1,Niters);
    curr_SDR_Linear = zeros(1,Niters);
    
    currSNR = snrLin(i);
    parfor n=1:Niters
        
        %% Generate source + side information
        x = sqrt(P)*randn(1,N);
        n_tilde = sqrt(P - rho^2)*randn(1,N);
        y = rho*x + n_tilde;
        w = (1/sqrt(currSNR))*randn(size(x));
        
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
        
        currenPower = mean(TxSig.^2);
        
        %% Decoding
        %% linear
        x_hat = ((rho/currSNR)*y + (1-rho^2)*RxSig_linear)/(1 + 1/currSNR - rho^2);
        curr_SDR_Linear(n) = sum((x_hat - x).^2);
        
        %% Tuncel Full MMSE
        sigmaW = sqrt(1/currSNR);
        x_hat_Tuncel = FullMMSE_decoder(RxSig_tuncel,y,alpha,beta,rho,sigmaW,codebook,centers,ones(size(probs)));
        
        curr_SDR_Tuncel_FullMMSE(n) = sum((x_hat_Tuncel - x).^2);
        
    end
    SDR_Linear_Final = sum(curr_SDR_Linear)/(N*Niters);
    SDR_Linear(i) = 10*log10(1/SDR_Linear_Final);
    
    SDR_Tuncel_FullMMSE_Final = sum(curr_SDR_Tuncel_FullMMSE)/(N*Niters);
    SDR_Tuncel_FullMMSE(i) = 10*log10(1/SDR_Tuncel_FullMMSE_Final);
    
    display(strcat('Finished SNR = ',num2str(snr(i))));
    
end
figure;hold all
plot(snr , 10*log10((1 + snrLin)/(1 - rho^2)),'LineWidth',2)
plot(snr , SDR_Tuncel_FullMMSE,'-kd','LineWidth',2)
plot(snr , 10*log10((1 + (1-rho^2)*snrLin)/(1-rho^2)),'r--','LineWidth',2)
% plot(snr , 10*log10((1 + (1-rho^2)*beta^2*snrLin)/(1-rho^2)),'--','LineWidth',2)

grid on; grid minor;
xlabel('SNR [dB]'); ylabel('SDR [dB]');
legend('SDR OPTA - SI case : (1 + SNR)/(1-\rho^2)',...
    strcat('Tuncel (MMSE Decoder) : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha),' \beta = ',num2str(beta)),...
    'Linear Scheme');
title(strcat('Distortion Vs SNR, \rho = ',num2str(rho)));

figure;hold all
plot(1./((1 + snrLin)/(1 - rho^2)),0.5*log(1 + snrLin),'LineWidth',2)
plot(1./(10.^(SDR_Tuncel_FullMMSE/10)),0.5*log(1 + snrLin),'-kd','LineWidth',2)
plot(1./((1 + (1-rho^2)*snrLin)/(1-rho^2)),0.5*log(1 + snrLin),'r--','LineWidth',2)
% plot(snr , 10*log10((1 + (1-rho^2)*beta^2*snrLin)/(1-rho^2)),'--','LineWidth',2)

grid on; grid minor;
xlabel('D '); ylabel('R');
legend('SDR OPTA - SI case : (1 + SNR)/(1-\rho^2)',...
    strcat('Tuncel (MMSE Decoder) : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha),' \beta = ',num2str(beta)),...
    'Linear Scheme');
title(strcat('Distortion Vs SNR, \rho = ',num2str(rho)));

figure;hold all
plot(10*log10(1./((1 + snrLin)/(1 - rho^2))),0.5*log(1 + snrLin),'LineWidth',2)
plot(10*log10(1./(10.^(SDR_Tuncel_FullMMSE/10))),0.5*log(1 + snrLin),'-kd','LineWidth',2)
plot(10*log10(1./((1 + (1-rho^2)*snrLin)/(1-rho^2))),0.5*log(1 + snrLin),'r--','LineWidth',2)
% plot(snr , 10*log10((1 + (1-rho^2)*beta^2*snrLin)/(1-rho^2)),'--','LineWidth',2)

grid on; grid minor;
xlabel('D[dB]'); ylabel('R');
legend('SDR OPTA - SI case : (1 + SNR)/(1-\rho^2)',...
    strcat('Tuncel (MMSE Decoder) : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha),' \beta = ',num2str(beta)),...
    'Linear Scheme');
title(strcat('Distortion Vs SNR, \rho = ',num2str(rho)));

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
x = (-7:dx:7);

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
