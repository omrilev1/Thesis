clear all; clc;
% test some code on the BEC

p = 0.35:0.01:0.45;

% generate code
[H,sortIndices] = parse_alist('.\Gallager\rate_0.5\MacKay_20000.r0.5','alist',0);% parse_alist('.\Irregular_optimized\BEC\rate_0.5\tenBrink__16383.v16.c8','alist',0);
ldpc_dec = comm.LDPCDecoder(H,'DecisionMethod','Soft decision','MaximumIterationCount',63,...
    'OutputValue','Whole codeword','NumIterationsOutputPort',1,'IterationTerminationCondition','Parity check satisfied');
nIter = 1e2;


channel = 'BEC' ; % 'BSC','BEC'
BER = zeros(size(p));
BLER = zeros(size(p));
for i=1:length(p)
    
    currBER = 0;
    currBLER = 0;
    for k=1:nIter
        
        ci = zeros(1,size(H,2));
        
        channelPattern = (rand(1,size(H,2)) > (1-p(i)));
        
        if strcmp(channel,'BEC')
            LLR = 100*ones(size(channelPattern));
            LLR(channelPattern == 1) = 0;
        else
            LLR = -1*log((1-p(i))/p(i))*ones(size(channelPattern));
            LLR(channelPattern == 0) = log((1-p(i))/p(i));
        end
        
        % decode LDPC
        [ldpc_decodedLLRs,numIterations] = step(ldpc_dec,LLR(:));
        
        % BER is the whole bits at LDPC decoder output with LLR = 0
        
        % count BER
        currBER = currBER + sum(ldpc_decodedLLRs == 0);
        currBLER = currBLER + (sum(ldpc_decodedLLRs(ldpc_decodedLLRs == 0)) > 0);
        
        
    end
    BER(i) = currBER / nIter / size(H,2);
    BLER(i) = currBLER / nIter;
end
BER(BER == 0) = 0.5/(nIter * size(H,2));
BLER(BLER == 0) = 0.5/(nIter);

figure;
semilogy(p,BER,'-o','LineWidth',2);
grid on; grid minor;
xlabel('Pr(Erasure)'); ylabel('Pr');
legend('BER');
title('ragular code, rate = 0.5 n = 20000')



