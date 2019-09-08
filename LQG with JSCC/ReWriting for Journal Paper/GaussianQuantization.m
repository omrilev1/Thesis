% Testing several transmission schemes of gaussian source where the
% receiver is equipped with side information
% - Regular Scaling and analog transmission
% - Tuncel Scheme, where we send the quantized value + quantization error,
% and use MMSE decoder in the reciever
% - Tuncel scheme, where the receiver knows T in a "genie aided" manner.
% (Upper bound on the optimal performance)
% - Modified Tuncel scheme, where we send T in a perfect "bit pipe
% channel", and send S in regular awgn. The overall rate is still the same.
% Note that for this scheme to work in "zero delay" we need to quantize for
% fixed length, whereas "fractional number of bits quantization" will
% result for a need to entropy encoder, and long delay

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
        codebook = generateGaussianCodeBook(partition);
end
nbits = calcBitsEntropy(3,delta);

N = 1;
P = 1;
rho = 0.9;
snr = (7:1:16);
alpha = 0.8;% 0.77;
beta = -1.15;% -1.1;

delta_kochman = 2.25;
beta_kochman = 1.05;

centers = [partition(1) (partition(1:end-1) + partition(2:end))/2 partition(end)];
Niters = 1e4;

SDR_Tuncel = zeros(size(snr));
SDR_Kochman = zeros(size(snr));
SDR_Linear = zeros(size(snr));
SDR_TuncelGenie = zeros(size(snr));
SDR_Hybrid_TuncelBitPipe = zeros(size(snr));

Bias_Tuncel = zeros(size(snr));
Rho_Tuncel = zeros(size(snr));

snrLin = 10.^(snr/10);
for i=1:length(snr)
    curr_SDR_Tuncel = 0;
    curr_SDR_Linear = 0;
    curr_SDR_TuncelGenie = 0;
    curr_SDR_Kochman = 0;
    curr_SDR_Hybrid_TuncelBitPipe = 0;
    
    curr_Bias_Tuncel = 0;
    curr_rho_Tuncel = 0;
    Tdist = zeros(1,Niters); Sdist = zeros(1,Niters);
    
    % calculate the SNR on S part, where we sacrifice rate for the integer
    % part
    hybridTuncelSNR = solveSNR(nbits,snrLin(i));
    errTuncel = zeros(1,Niters);
    errGenie = zeros(1,Niters);
    TuncelK = zeros(1,Niters);
    for n=1:Niters
        
        % source + side information
        x = sqrt(P)*randn(1,N);
        y = rho*x + sqrt(P - rho^2)*randn(1,N);
        
        %% Encoding + Channel
        % linear
        RxSig_linear = x + (1/sqrt(snrLin(i)))*randn(size(x));
        
        % Tuncel
        currDistance = abs(x - centers(:));
        [~,minIdx] = min(currDistance,[],1);
        T = codebook(minIdx);
        S = x - T;
        
        TxSig = alpha*T + beta*S;
        RxSig_tuncel = TxSig + sqrt(1/snrLin(i))*randn(size(TxSig));
        
        % Debug - Manual Power Calibration
        %                                 currenPower = mean(TxSig.^2);
        %                                 figure;
        %                                 subplot(221); histogram(S); title(strcat('S part histogram, E(s^2) = ',num2str(mean(S.^2))));grid on; grid minor;
        %                                 subplot(222); histogram(TxSig); title(strcat('TxSig histogram, E(TxSig^2) = ',num2str(mean(TxSig.^2))));grid on; grid minor;
        %                                 subplot(223); scatter(rho*x(1:10:end),TxSig(1:10:end)); title('[TxSig,rho*x]');xlabel('rho*x');ylabel('TxSig');grid on; grid minor;
        %                                 subplot(224); scatter(y(1:10:end),RxSig_tuncel(1:10:end)); title(strcat('[RxSig,Y], SNR = ',num2str(snr(i)),'[dB]'));xlabel('y');ylabel('RxSig');grid on; grid minor;
        %
        
        % Kochman Zamir
        curr_dither = 0; % 2*delta_kochman*rand - delta_kochman;
        TxSig = mod(beta_kochman*x + delta_kochman + curr_dither,2*delta_kochman) - delta_kochman;
        RxSig_kochman = TxSig + sqrt(1/snrLin(i))*randn(size(TxSig));
        
%         % Debug - Manual Power Calibration
%               currenPower = mean(TxSig.^2);
%         figure;
%         subplot(211); histogram(TxSig); title(strcat('TxSig histogram, E(TxSig^2) = ',num2str(mean(TxSig.^2))));grid on; grid minor;
%         subplot(212); scatter(rho*x(1:10:end),TxSig(1:10:end)); title('[TxSig,rho*x]');xlabel('rho*x');ylabel('TxSig');grid on; grid minor;
%         subplot(213); scatter(y(1:10:end),RxSig_kochman(1:10:end)); title(strcat('[RxSig,Y], SNR = ',num2str(snr(i)),'[dB]'));xlabel('y');ylabel('RxSig');grid on; grid minor;
%         
%         
%   
        
        % Hybrid Tuncel BitPipe (Infinite Delay)
        
        [currPower,~] = calcQuantizationErrorPower(T,codebook,delta);
        scale = 1/sqrt(currPower);
        RxSig_tuncel_BitPipe = sqrt(scale)*S + sqrt(1/hybridTuncelSNR)*randn(size(TxSig));
        
        %% Decoding
        % linear
        x_hat = ((rho/snrLin(i))*y + (1-rho^2)*RxSig_linear)/(1 + 1/snrLin(i) - rho^2);
        curr_SDR_Linear = curr_SDR_Linear + sum((x_hat - x).^2);
        
        % Tuncel
        sigmaW = sqrt(1/snrLin(i));
        T_hat = MMSE_decoder_Int(RxSig_tuncel,y,delta,codebook,alpha,beta,rho,sigmaW,NaN);
        
        S_hat = MMSE_decoder_Frac(RxSig_tuncel,y,delta,alpha,beta,rho,sigmaW,codebook,NaN);
        x_hat_Tuncel = T_hat + S_hat;
        
        curr_SDR_Tuncel = curr_SDR_Tuncel + sum((x_hat_Tuncel - x).^2);
        curr_Bias_Tuncel = curr_Bias_Tuncel + x_hat_Tuncel - x ;
        curr_rho_Tuncel = curr_rho_Tuncel + x*(x_hat_Tuncel - x);
        
        % Tuncel Genie - the receiver knows T
        S_hat = MMSE_decoder_Frac(RxSig_tuncel,y,delta,alpha,beta,rho,sigmaW,codebook,T);
        x_hat_Tuncel_Genie = T + S_hat;
        
        curr_SDR_TuncelGenie = curr_SDR_TuncelGenie + sum((x_hat_Tuncel_Genie - x).^2);
        
        % Kochman
        x_hat_kochman = Kochman_Decoder(RxSig_kochman,y,beta_kochman,delta_kochman,rho,sigmaW,curr_dither);
        
        curr_SDR_Kochman = curr_SDR_Kochman + sum((x_hat_kochman - x).^2);
        
        % Tuncel Hybrid, infinite delay
        sigmaW_hybrid = sqrt(1/hybridTuncelSNR);
        S_hat = MMSE_decoder_Frac(RxSig_tuncel_BitPipe + T,y,delta,1,sqrt(scale),rho,sigmaW_hybrid,codebook,T);
        x_hat_Tuncel_Hybrid = T + S_hat;
        
        curr_SDR_Hybrid_TuncelBitPipe = curr_SDR_Hybrid_TuncelBitPipe + sum((x_hat_Tuncel_Hybrid - x).^2);
        
        errTuncel(n) = x_hat_Tuncel - x;
        errGenie(n) = x_hat_Tuncel_Genie - x;
        TuncelK(n) = TxSig - 0.2580*x;
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
    
    SDR_TuncelBitPipe_Final = curr_SDR_Hybrid_TuncelBitPipe/(N*Niters);
    SDR_Hybrid_TuncelBitPipe(i) = 10*log10(1/SDR_TuncelBitPipe_Final);
    
    display(strcat('Finished SNR = ',num2str(snr(i))));
end
figure;hold all
plot(snr , 10*log10((1 + snrLin)/(1 - rho^2)),'LineWidth',2)
plot(snr , SDR_TuncelGenie,'g.-','LineWidth',2)
plot(snr , SDR_Tuncel,'-kd','LineWidth',2)
plot(snr , SDR_Kochman,'--gs','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
plot(snr , SDR_Linear,'r--','LineWidth',2)
plot(snr , SDR_Hybrid_TuncelBitPipe,'b-p','LineWidth',2)
grid on; grid minor;
xlabel('SNR [dB]'); ylabel('SDR [dB]');
legend('SDR OPTA - SI case : (1 + SNR)/(1-\rho^2)',...
    strcat('Tuncel (MMSE Decoder) With Genie : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha),' \beta = ',num2str(beta)),...
    strcat('Tuncel (MMSE Decoder) : \Delta = ',num2str(delta),' \alpha = ',num2str(alpha),' \beta = ',num2str(beta)),...
    strcat('Kochman - Zamir : \Delta = ',num2str(delta_kochman),' \beta = ',num2str(beta_kochman)),...
    'Linear Scheme',...
    strcat('Hybrid Tuncel + BitPipe : \Delta = ',num2str(delta),' n_{bits} = ',num2str(nbits)));
title(strcat('Distortion Vs SNR, \rho = ',num2str(rho)));
ylim([10 20]);


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

function [s_hat] = Kochman_Decoder(RxSig,y,beta,delta,rho,sigmaW,curr_dither)

mod_sig = mod(RxSig - beta*rho*y - curr_dither + delta,2*delta) - delta;

alpha_s = (beta*(1-rho^2)) / (beta^2 * (1-rho^2) + sigmaW^2);
s_hat = rho*y + alpha_s * mod_sig;

end

function [codebook] = generateGaussianCodeBook(partition)

dx = 0.001;
x = (-3 + partition(1)):dx:(partition(end) + 3);

pdf = (1/sqrt(2*pi))*exp(-0.5*x.^2);
codebook = zeros(1,length(partition) - 1);
for i=1:(length(partition) - 1)
    
    curIdx = find((x >= partition(i)) & (x <= partition(i+1)));
    curPDF = pdf(curIdx) / sum(pdf(curIdx));
    curX = x(curIdx);
    codebook(i) = sum(curX.*curPDF);
end

% generate first and last codewords

% first
curIdx = find(x <= partition(1));
curPDF = pdf(curIdx) / sum(pdf(curIdx));
curX = x(curIdx);
codebook = [sum(curX.*curPDF) codebook];

% last codeword
curIdx = find(x >= partition(end));
curPDF = pdf(curIdx) / sum(pdf(curIdx));
curX = x(curIdx);
codebook = [codebook sum(curX.*curPDF)];
end

function [hybridTuncelSNR] = solveSNR(int_capacity,snrLin)

effectiveCapacity = 0.5*log2(1+snrLin) - int_capacity;
if effectiveCapacity < 0
    hybridTuncelSNR = -inf;
else
    hybridTuncelSNR = 2^(2*effectiveCapacity) - 1;
end
end

function [nbits] = calcBitsEntropy(nSegments,delta)

% calculate the exact number of bits we need in order to send the quantized
% source, after quantizing to nSegments bits
dx = 0.001;
x = (-nSegments*delta/2):dx:(nSegments*delta/2);

pdf = (1/sqrt(2*pi))*exp(-0.5*x.^2);
probs = zeros(1,nSegments);
for i=1:nSegments
    
    if i == 1
        probs(i) = sum(pdf((x <= -(delta/2)*(nSegments-1)/2)))/ sum(pdf);
    elseif i==nSegments
        probs(i) = sum(pdf((x >= (delta/2)*(nSegments-1)/2)))/ sum(pdf);
    else
        probs(i) = sum(pdf((x >= -(delta/2)*(nSegments-1)/2 + (i-2)*delta) &...
            (x <= -(delta/2)*(nSegments-1)/2 + (i-1)*delta)))/ sum(pdf);
    end
end
logProbs = log2(probs);
logProbs(logProbs == -inf) = 0;
nbits = -1*sum(probs.*logProbs);
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
