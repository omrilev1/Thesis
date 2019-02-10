clear all; close all; clc;
addpath(genpath('.\TimeVarying_ConvCode'))
addpath(genpath('.\QAM'))
addpath(genpath('.\results'))

% This script simulate concatenated coding scheme, where the outer code is
% high rate LDPC, and the inner code is time varying convolutional code
% Every thing implemented using matlab legacy functions
% We compare to :
% LTE Turbo
% 5G LDPC

resultsFileName = strcat('.\results\',date,'DVBS2_LDPC');


% Init Params

% the next rates seems to not work : [1/2 1/3 3/5 8/9]
dvb_rate = [1/4 2/5 2/3 3/4 5/6 9/10]; 
cnt = 1;

EbN0 = (-0.5:0.25:7);
m = 1;

Nout = 50;
Nin = 50;

BER = zeros(length(EbN0),size(dvb_rate,2));

NbitsLDPC = 64800;
while cnt<=size(dvb_rate,2)
    curr_LDPC_r = dvb_rate(cnt);
    
    % generate LDPC H matrix
    H = dvbs2ldpc(curr_LDPC_r);
    ldpc_enc = comm.LDPCEncoder(H);
    ldpc_dec = comm.LDPCDecoder(H,'DecisionMethod','Soft decision','MaximumIterationCount',70);
    
    Nbits = size(H,2) - size(H,1);numOfErr = 0;
    for i=1:length(EbN0)
        numOfErr = 0;
        for n_in = 1:Nin
            SNR = EbN0(i) + 10*log10(curr_LDPC_r) + 10*log10(m);
            raw_bits = rand(Nbits,1) > 0.5;
            
            %% Encode LDPC
            Tx_codeword = step(ldpc_enc,raw_bits);
            
            %% channel
            x = qam_modulate(Tx_codeword(:),m,1);
            
            if m==1
                x = x.*exp(-1j*pi/4);
            end
            y = x + (1/sqrt(2))*10^(-SNR/20)*(randn(size(x)) + 1j*randn(size(x)));
            
            llr = calcLLR(y,m,10^(-SNR/10));
            
            %% decode LDPC
            ldpc_decodedLLRs = step(ldpc_dec,llr);
            ldpc_decodedBits = ldpc_decodedLLRs > 0;
            
            % count BER
            numOfErr = numOfErr + sum(ldpc_decodedBits(:) ~= raw_bits);
            
            if numOfErr > 50
                BER(i,cnt) = numOfErr/(Nbits*n_in);
                break
            end
        end
        if BER(i,cnt) < 0.5*1e-6
            BER(i:end,cnt) = 0.5*1e-6;
            
            iterMessage = sprintf('LDPC : rate %0.2f , EbN0 = %2g',curr_LDPC_r,EbN0(i));
            iterMessage2 = sprintf('BER = %.5f',BER(i,cnt));
            disp(iterMessage);disp(iterMessage2);disp('%%%%%%%%%%%%%%%%%%%%%');
            break
        end
        
        iterMessage = sprintf('LDPC : rate %0.2f , EbN0 = %2g',curr_LDPC_r,EbN0(i));
        iterMessage2 = sprintf('BER = %.5f',BER(i,cnt));
        disp(iterMessage);disp(iterMessage2);disp('%%%%%%%%%%%%%%%%%%%%%');
    end
    cnt = cnt + 1;
end


lineType = ['-kd';'-go';'-b*';'-rs';'--c';'-m+'];
figure;
for idx = 1:(length(dvb_rate))
    semilogy(EbN0.',BER(:,idx),lineType(idx,:),'LineWidth',2)
    hold on;
end
grid on; grid minor;
xlabel('E_b / N_0 [dB]'); ylabel('BER')

legend(sprintf('Rate = %0.2f',dvb_rate(1)),...
    sprintf('Rate = %0.2f',dvb_rate(2)),...
    sprintf('Rate = %0.2f',dvb_rate(3)),...
    sprintf('Rate = %0.2f',dvb_rate(4)),...
    sprintf('Rate = %0.2f',dvb_rate(5)),...
    sprintf('Rate = %0.2f',dvb_rate(6)));

plotMessage = 'Plain Vanilla DVB-S2 LDPC';
title(plotMessage)

save(resultsFileName,'dvb_rate','BER');

