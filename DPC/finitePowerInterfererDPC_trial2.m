%% Evaluate the MI in various DPC problems
close all; clear; clc;

addpath(genpath('Utils'));



SNR = -10:2:10;
SIR = -3;
SNRlin = 10.^(SNR/10);
SIRlin = 10.^(SIR/10);

% generate histogram with ~sqrt(N) bins
frmLen = 1e5; Niter = 1e3;

% Init Arrays
SER = zeros(length(SNR),length(SIR));

inputMod = 'uni'; % 'uni' for uniform signal , 'mPAM' for m order PAM signal
pamOrder = log2(2);
intPdf = 'gaussian';

% define parameters
Delta = (2^pamOrder-1) + 2; % lattice constant , 1D

% Constellation Prior Distribution : calculated empirically in advanced
switch pamOrder
    case 1
        % Estimated Prior : obtained from simulation and using laplace correction
        Prior = [5.5*1e-7 3.6645*1e-4 0.0261 0.2438 0.4594 0.2438 0.0261 3.6515*1e-4 7.5*1e-7];
        Prior = Prior/sum(Prior);
        % Segment relevant to base segment : segment i means the original
        % symbol obtained by subtracting i*Delta
        Segment = [-4 -3 -2 -1 0 1 2 3 4];
        
        % generate table of original symbols (of augmented constellation)
        % and corresponding probabilities
        origSymbols = [-1 1];
        augmentedSymbols = [];
        augmentedSymbolsPrior = [];
        for i=1:length(Segment)
            augmentedSymbols = [augmentedSymbols  origSymbols + Segment(i)*2*Delta];
            pamPrior = 1/(2*pamOrder);
            augmentedSymbolsPrior = [augmentedSymbolsPrior pamPrior*ones(1,2)*Prior(i)];
        end
end

% Here Statistics loop starts
for i=1:length(SNR)
    for j=1:length(SIR)
        numOfErrs = zeros(1,Niter);
        for iterIn = 1:Niter
            
            % generate random stream  : We start from regular PAM
            ai = bi2de(reshape((rand(1,frmLen*pamOrder) > 0.5),[],pamOrder));
            pamSymbols = pammod(ai,2^pamOrder);
            in = pamSymbols; % pamSymbols + 2^pamOrder - 1 + 0.05;
            % calculate PAM index in constellation , just for decoding
            inIdx = ones(size(in));
            inIdx(in > 0) = 2;
            
            % calculate powers to calibrate simulation SNR
            Px = mean(in.^2);
            sigmaS = sqrt(Px/SIRlin(j));
            
            % generate Interferer and calculate Tx Signal
            switch intPdf
                case 'gaussian'
                    S = sigmaS*randn(frmLen,1);
                case 'PAM'
                    ai = bi2de(reshape((rand(1,frmLen*pamOrder) > 0.5),[],pamOrder));
                    interfererPamSymbols = pammod(ai,2^pamOrder);
                    S = interfererPamSymbols + 2^pamOrder - 1;
                    S = sigmaS * S/sqrt(mean(S.^2));
            end
            
            TxOutput = mod(in - S + Delta,2*Delta) - Delta;
            sigmaZ = sqrt(mean(TxOutput.^2)/SNRlin(i));
            
            % generate wgn noise
            Z = sigmaZ*randn(frmLen,1);
            
            % generate channel output
            RxInput = TxOutput + Z + S;
            
            % Decode Original PAM symbol and calculate number of errors :
            % squared error between augmented constellation and recieved
            % sample + log(prior)
            
            logLikelihoodFunction = (1/sigmaZ) * abs(RxInput(:) - augmentedSymbols).^2 ...
                - log(repmat(augmentedSymbolsPrior,length(RxInput),1));
            
            [~,decodedSymbolIdx] = min(logLikelihoodFunction,[],2);
            decodedSymbolIdx = mod(decodedSymbolIdx,2);
            decodedSymbolIdx(decodedSymbolIdx == 0) = 2;
            
            % Count Errors
            numOfErrs(iterIn) = sum(decodedSymbolIdx ~= inIdx);
            
            % break if achieve enough statistics
            if(sum(numOfErrs) > 10)
                break
            end
        end
        SER(i,j) =  sum(numOfErrs)/iterIn/frmLen;
    end
    if (SER(i,j) < 1e-6)
        SER(SER(:,j) == 0) = 0.5/(frmLen * Niter);
        break
    end
end
%% Plot Capacity(SNR) and Capacity(Eb/N0)
PAM2SER = qfunc(sqrt(2*10.^((SNR - 5)/10)));

% Capacity(SNR)
figure;
semilogy(SNR,PAM2SER,'--','LineWidth',1.5);hold on;
semilogy(SNR,SER,'-*','LineWidth',1.5)
grid on; grid minor;
xlabel('SNR [dB]'); ylabel('Probability Of Symbol Error');
title('2PAM With Interferer Symbol Error Probability')

