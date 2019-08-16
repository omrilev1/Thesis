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
dvb_rate = [1/4 1/3 1/2 2/5 2/3 3/4 5/6 8/9];
cnt = 1;

m = 1;

cap = 0.35:0.03:0.47;

Nout = 10;
Nin = 5;

type = 'long';
rates = [2/5]; % 1/4 , 2/3 

channelType = [2,1,0]; % {0 == 'BEC',1 == 'BSC',2 == 'AWGN'};

NbitsLDPC = 64800;
params = combvec(channelType,rates);

BER = zeros(length(cap),size(params,2));
while cnt<=size(params,2)
    
    curr_LDPC_r = params(2,cnt);
    curr_chan = params(1,cnt);
    % generate LDPC H matrix
    H = dvbs2ldpc(curr_LDPC_r);
    
    ldpc_enc = comm.LDPCEncoder(H);
    ldpc_dec = comm.LDPCDecoder(H,'DecisionMethod','Soft decision','MaximumIterationCount',63,'IterationTerminationCondition','Parity check satisfied');
    
    
    Nbits = size(H,2) - size(H,1);numOfErr = 0;
    for i=1:length(cap)
        numOfErr = 0;
        for n_out = 1:Nout
            for n_in = 1:Nin
                
                raw_bits = rand(Nbits,1) > 0.5;
                
                %% Encode LDPC
                Tx_codeword = step(ldpc_enc,raw_bits);
                
                %% channel
                x = qam_modulate(Tx_codeword(:),m,1);
                
                if m==1
                    x = real(x.*exp(-1j*pi/4));
                end
                
                [llr] = channelModel(x,cap(i),curr_chan);
                
                
                %% decode LDPC
                ldpc_decodedLLRs = step(ldpc_dec,llr);
                ldpc_decodedBits = ldpc_decodedLLRs > 0;
                
                % count BER
                numOfErr = numOfErr + sum(ldpc_decodedBits(:) ~= raw_bits);
                
            end
            
            if numOfErr > 50
                BER(i,cnt) = numOfErr/(Nbits*n_in);
                break
            end
        end

        if BER(i,cnt) == 0
            BER(i,cnt) = 0.5/(Nbits*Nin*Nout);
        end
        iterMessage = sprintf('LDPC : rate %0.2f , capacity = %2g',curr_LDPC_r,cap(i));
        iterMessage2 = sprintf('BER = %.5f',BER(i,cnt));
        disp(iterMessage);disp(iterMessage2);disp('%%%%%%%%%%%%%%%%%%%%%');
    end
    cnt = cnt + 1;
    if mod(cnt-1,3) == 0
        resultPlot(cap,BER,3,curr_LDPC_r,cnt)
    end
    
end

save(resultsFileName,'dvb_rate','BER');

function [] = resultPlot(capacity,BER,numOfLines,rLDPC,cnt)
% plot group of channels
lineType = ['-kd';'-go';'-b*';'-rs';'--c';'-xm'];
figure;
for idx = 1:numOfLines
    semilogy(capacity.',BER(:,numOfLines*((cnt - 1)/numOfLines - 1) + idx),lineType(mod(idx,length(lineType)),:),'LineWidth',2)
    hold on;
end
semilogy(rLDPC*ones(1,100),linspace(1e-4,1,100),'--k','LineWidth',1.5)
grid on; grid minor;
xlabel('C [nats]'); ylabel('BER')
legend('AWGN','BSC','BEC')
title(strcat(' R_{LDPC} = ',num2str(rLDPC)))
end


