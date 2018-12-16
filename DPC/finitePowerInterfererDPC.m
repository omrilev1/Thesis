%% Evaluate the MI in various DPC problems
% This script calculate equivalent noise pdf and MI , for various of input
% distributions and interferer distributions

close all; clear; clc;
isHistogramPlotRun = 1;

addpath(genpath('Utils'));
% DPC with finitw power interferer : we only take modulo at transmitter
calcMethod = 'nnEntropy';


SNR = 100 + (10:2:20); % -5:2:10;
SIR = 0;
SNRlin = 10.^(SNR/10);
SIRlin = 10.^(SIR/10);


% generate histogram with ~sqrt(N) bins
NforStatistics = 1e7;
NforMI = round(sqrt(NforStatistics));
% Init Arrays
IxyAWGN = zeros(length(SNR),length(SIR));

% input modulations and interferer type
inputMod = {'PAM','uni'}; % 'uni' for uniform signal , 'mPAM' for m order PAM signal
intPdf = {'gaussian','uni','PAM'};
pamOrder = 2;

% define parameters
Delta = pamOrder - 0.001; % lattice constant , 1D
Px = 1;

for i=1:length(SNR)
    for j=1:length(SIR)
        
        %% generate input signal  : uniform between [-Delta,Delta] and M-PAM : [-(M-1),-(M-3),...,-1,1,...,M-3,M-1]
        % uniform input : between -delta to delta
        uniformInput = Delta * 2 * rand(NforStatistics,1) - Delta;
        
        % pam input
        ai = bi2de(reshape((rand(1,NforStatistics*log2(pamOrder)) > 0.5),[],log2(pamOrder)));
        pamInput = pammod(ai,pamOrder);
        
        InputMat = [uniformInput(:) pamInput(:) pamInput(:) pamInput(:)];
        
        %% generate Interferer
        % gaussian interferer with the desired SIR
        sigmaS = sqrt(mean(pamInput.^2)/SIRlin);
        gaussianInterferer =  sigmaS*randn(NforStatistics,1);
        
        % PAM interferer
        ai = bi2de(reshape((rand(1,NforStatistics*log2(pamOrder)) > 0.5),[],log2(pamOrder)));
        pamInterferer = pammod(ai,pamOrder);
        
        % uniform interderer
        uniformInterferer = Delta * 2 * rand(NforStatistics,1) - Delta;
        
        InterfererMat = [uniformInterferer(:) uniformInterferer(:) gaussianInterferer(:) pamInterferer(:)];
        
        %% generate Tx : subtract interferer and modulo
        
        % uniform input , uniform interferer
        TxOutput_UniUni = mod(uniformInput - uniformInterferer + Delta,2*Delta) - Delta;
        
        % pam input , uniform interferer
        TxOutput_PamUni = mod(pamInput - uniformInterferer + Delta,2*Delta) - Delta;
        
        % pam input , gaussian interferer
        TxOutput_PamGauss = mod(pamInput - gaussianInterferer + Delta,2*Delta) - Delta;
        
        % pam input , pam interferer
        TxOutput_PamPam = mod(pamInput - pamInterferer + Delta,2*Delta) - Delta;
        
        TxSigMat = [TxOutput_UniUni(:) TxOutput_PamUni(:) TxOutput_PamGauss(:) TxOutput_PamPam(:)];
        
        %% Channel : noise + interferer
        % calibrate the noise variance ccording to Tx Output signal - it
        % should be almost the same for all kind of modulations. We
        % calibrate accroding to the worst case (higher power)
        Px = max(mean(TxSigMat.^2,1));
        sigmaZ = sqrt(Px/SNRlin(i));
        
        % generate wgn noise
        Z = sigmaZ*randn(size(TxSigMat));
        
        % generate channel output
        RxInput = TxSigMat + Z + InterfererMat;
        
        %% Calculate statistics : pdf's of error and MI
        % Calculate entropy and MI
        errSig = RxInput - InputMat;
        augmentedConstellation = TxSigMat + InterfererMat;
        
        if isHistogramPlotRun
            % pdf of error signal
            figure;
            [pdf_UniUni] = histogram(errSig(:,1),'Normalization','probability','NumBins',2*round(sqrt(NforMI)));
            grid on; grid minor; title('Pdf of Error signal , Uniform Input and Uniform Interferer');
            legend('SNR  = 10[dB]')

            figure;
            [pdf_PamUni] = histogram(errSig(:,2),'Normalization','probability','NumBins',2*round(sqrt(NforMI)));
            grid on; grid minor; title('Pdf of Error signal , PAM Input and Uniform Interferer');
            legend('SNR  = 10[dB]')
            
            figure;
            [pdf_PamGauss] = histogram(errSig(:,3),'Normalization','probability','NumBins',2*round(sqrt(NforMI)));
            grid on; grid minor; title('Pdf of Error signal , PAM Input and Gaussian Interferer');
            legend('SNR  = 10[dB]')
            
            figure;
            [pdf_PamPam] = histogram(errSig(:,4),'Normalization','probability','NumBins',2*round(sqrt(NforMI)));
            grid on; grid minor; title('Pdf of Error signal , PAM Input and PAM Interferer');
            legend('SNR  = 10[dB]')
            
            % pdf of augmented constellation
            figure;
            [inputPDF_UniUni] = histogram(augmentedConstellation(:,1),'Normalization','probability','NumBins',round(sqrt(NforMI)));
            grid on; grid minor; title('Pdf of Augmented Symbol , Uniform Input and Uniform Interferer');
            legend('SIR  = 0[dB]')
            
            figure;
            [inputPDF_PamUni] = histogram(augmentedConstellation(:,2),'Normalization','probability','NumBins',100);
            grid on; grid minor; title('Pdf of Augmented Symbol , PAM Input and Uniform Interferer');
            legend('SIR  = 0[dB]')
            
            figure;
            [inputPDF_PamGauss] = histogram(augmentedConstellation(:,3),'Normalization','probability','NumBins',2*round(sqrt(NforMI)));
            grid on; grid minor; title('Pdf of Augmented Symbol , PAM Input and Gaussian Interferer');
            legend('SIR  = 0[dB]')
            
            figure;
            [inputPDF_PamPam] = histogram(augmentedConstellation(:,4),'Normalization','probability','NumBins',10);
            grid on; grid minor; title('Pdf of Augmented Symbol , PAM Input and PAM Interferer');
            legend('SIR  = 0[dB]')
        else
            [Hx,Hy,Hxy,IxyAWGN(i,j)] = mutualInfoCalculator(in,RxInput,NforMI,calcMethod);
        end
    end
end

%% Plot Capacity(SNR) and Capacity(Eb/N0)
awgnCapTheory = 0.5*log2(1 + SNRlin);

% Capacity(SNR)
figure;
plot(SNR,awgnCapTheory,'--','LineWidth',1.5);hold on;
plot(SNR,IxyAWGN,'-*','LineWidth',1.5);hold on;
plot(SNR,IxyAWGNuni,'-o','LineWidth',1.5);hold on;
plot(SNR,IxyDPC,'-^','LineWidth',1.5);hold on;
grid on; grid minor;
xlabel('SNR [dB]'); ylabel('I(X;Y)');
title('Capacity in various scenarios')
legend('Theoretic AWGN','Simulation AWGN','AWGN with Uniform input','DPC');

% Capacity(Eb/N0)
figure;
plot(10*log10(SNRlin./(2*awgnCapTheory)),awgnCapTheory,'--','LineWidth',1.5);hold on;
plot(10*log10(SNRlin./(2*IxyAWGN)),IxyAWGN,'-*','LineWidth',1.5);hold on;
grid on; grid minor;
xlabel('SNR [dB]'); ylabel('I(X;Y)');
title('Capacity in various scenarios')
legend('Theoretic AWGN','AWGN With Uniform Input','DPC Without Modulo','DPC With Modulo');
