clear all; clc;
addpath(genpath('.\QAM'))
addpath(genpath('.\results'))
addpath(genpath('.\Ensembles'))

%% This script simulate concatenated coding scheme :
% The outer code is high rate LDPC. The inner code suppose to be
% some code with performance guarantee over BSC (for example convolutional codes or polar codes).
% The main thing we check is the universality of the code

% The overall coding scheme is linear, thus we simulate only the zeros codeword
% and randomize the noise

% The simulation can be done for AWGN, BSC and BEC

resultsFileName = strcat('.\results\',date,'GallagerEnsemble');

% Init Params
convCode.constraintLength = 10;
convCode.rate = 1/2;
convCode.rsc = 1;

fullBP_iterations = 0;

decoding_method = 1; % 0 for LDPC in hard decision , 1 for soft

% LDPC Ensemble
Ensemble = 'Irregular_BEC';% 'Irregular_BEC';% 'Gallager';
ldpc_rate = 0.66;% [0.87,0.9,0.82,0.66];

puncture_percent = 0.19;% [0.22 0.15 0.1];
channelType = [2,1,0]; % {0 == 'BEC',1 == 'BSC',2 == 'AWGN'};

blockIntrlvDepth = 256;
params = combvec(fullBP_iterations,puncture_percent,decoding_method,channelType,ldpc_rate);

cnt = 1;

cap = 0.37:0.015:0.6;
m = 1;

Nout = 5;
Nin = 10;

BER = zeros(length(cap),size(params,2));
BER_conv = zeros(length(cap),size(params,2));
effectiveRate = zeros(1,size(params,2));

while cnt<=size(params,2)
    curr_iter = params(1,cnt);
    curr_punct = params(2,cnt);
    curr_decMethod = params(3,cnt);
    curr_ldpc_rate = params(5,cnt);
    curr_chan = params(4,cnt);
    %% Generate the corresponding convolutional code
    [poly,conv_enc,conv_dec] = convCodeGen(convCode.constraintLength,convCode.rate,convCode.rsc);
    
    % generate LDPC H matrix
    [H,curr_LDPC_r,ldpc_dec] = ldpcGen(curr_ldpc_rate,Ensemble,curr_decMethod);
    
    Nbits = size(H,2);numOfErr_Conv = 0; numOfErr = 0;
    Ndata = Nbits - size(H,1);
    for i=1:length(cap)
        numOfErr_Conv = zeros(1,curr_iter+1);
        numOfErr_Final = zeros(1,curr_iter+1);
        
        for n_out = 1:Nout
            for n_in = 1:Nin
                
                %% Encode LDPC
                ldpc_codedBits = zeros(1,Nbits);
                
                %% Add deletions and interleave bits
                numOfPunct = Nbits*curr_punct;
                effectiveConvLen = Nbits - numOfPunct;
                
                % truncate to use block interleaver - we add little more
                % deletions
                NbitsConv = blockIntrlvDepth*floor(effectiveConvLen/blockIntrlvDepth);
                numOfPunct =  numOfPunct + effectiveConvLen - NbitsConv ;
                
                % remove bits, that will be considered as deleted in the LDPC
                % decoding stage
                punctIdx = randperm(Nbits,numOfPunct);
                dataIdx = setdiff(1:Nbits,punctIdx);
                
                ldpc_codedBits_ForConv = ldpc_codedBits(dataIdx);
                ldpc_codedBits_ForConv = (ldpc_codedBits_ForConv > 0.5);
                
                blockIntrlvMat = reshape(ldpc_codedBits_ForConv,[],blockIntrlvDepth);
                blockIntrlv    = blockIntrlvMat.';
                interleavedBits = blockIntrlv(:);
                
                %% Encode convolutional
                % simple scrambler
                scramblingSeq = rand(size(interleavedBits)) > 0.5;
                scrambledBits = mod(scramblingSeq + interleavedBits,2);
                Tx_codeword = step(conv_enc,scrambledBits);
                
                x = qam_modulate(Tx_codeword(:),m,1);
                
                if m==1
                    x = real(x.*exp(-1j*pi/4));
                end
                
                %% channel
                
                currR = Ndata/length(Tx_codeword);
                [llr] = channelModel(x,cap(i),curr_chan);
                
                %% Decode Convolutional code + LDPC in an iterative manner
                
                convCode_priors = zeros(size(interleavedBits));
                for decodeIter=0:curr_iter
                    
                    % Decode Convolutional code
                    
                    % in the first iteraticon take aprior probabilities to be
                    % [1/2,1/2]. After that take the apriori probabilities from
                    % LDPC output LLRs
                    dataSoft = step(conv_dec,convCode_priors,llr);
                    dataHard = dataSoft > 0;
                    if curr_decMethod == 0
                        llr_ForDeintlrv = dataHard;
                    else
                        llr_ForDeintlrv = dataSoft;
                    end
                    llr_ForDeintlrv = llr_ForDeintlrv.*(2*scramblingSeq - 1);
                    
                    % count Convolutional code BER
                    numOfErr_Conv(decodeIter+1) = numOfErr_Conv(decodeIter+1) + sum(dataHard(2:end-1) ~= scrambledBits(2:end-1));
                    
                    %% Deinterleave soft bits
                    block_DeIntrlv = reshape(llr_ForDeintlrv,blockIntrlvDepth,[]);
                    block_DeIntrlvMat = block_DeIntrlv.';
                    LDPC_llrs = block_DeIntrlvMat(:);
                    
                    % add punctures - punctures will be treated as LLR = 0
                    LDPC_final_llrs = zeros(Nbits,1);
                    LDPC_final_llrs(dataIdx) = LDPC_llrs;
                    LDPC_final_llrs(punctIdx) = -0.01;
                    
                    % decode LDPC
                    [ldpc_decodedLLRs,numIterations] = step(ldpc_dec,LDPC_final_llrs);
                    ldpc_decodedBits = ldpc_decodedLLRs < 0;
                    
                    % count BER
                    if (sum(abs(ldpc_decodedLLRs - LDPC_final_llrs)) == 0)
                        numOfErr_Final(decodeIter+1) = numOfErr_Final(decodeIter+1) + Ndata/2;
                    else
                        numOfErr_Final(decodeIter+1) = numOfErr_Final(decodeIter+1) + sum(ldpc_decodedBits(2:Ndata-1));
                    end
                    
                    % calculate conv-code next iteration llrs
                    convCode_priors = ldpc_decodedLLRs(dataIdx);
                    blockIntrlvMat = reshape(convCode_priors,[],blockIntrlvDepth);
                    blockIntrlv    = blockIntrlvMat.';
                    convCode_priors = blockIntrlv(:);
                end
            end
            if numOfErr_Final(end) > 50
                BER(i,cnt) = numOfErr_Final(end)/((Ndata-2)*Nin);
                break
            end
        end
        
        BER_conv(i,cnt) = numOfErr_Conv(decodeIter+1)/(n_in*n_out*length(interleavedBits));
        effectiveRate(cnt) = Ndata/length(Tx_codeword);
        iterMessage1 = sprintf('LDPC : rate %0.2f , Punct = %0.2f \nConvolutional Rate = %0.2f constraint length = %0.2f \n Effective Rate = %0.2f%',...
            curr_LDPC_r,100*curr_punct,convCode.rate,convCode.constraintLength,effectiveRate(cnt));
        iterMessage2 = sprintf('BER Convolutional = %.5f overall BER = %.5f',BER_conv(i,cnt),BER(i,cnt));
        iterMessage3 = sprintf('C = %2g',cap(i));
        iterMessage4 = sprintf('Channel is %2g',curr_chan);
        disp(iterMessage3);disp(iterMessage4);disp(iterMessage1);disp(iterMessage2);disp('%%%%%%%%%%%%%%%%%%%%%');
    end
    
    cnt = cnt + 1;
    
    if mod(cnt-1,3) == 0
        resultPlot(cap,BER,3,length(Tx_codeword),curr_ldpc_rate,0.5,Ndata/length(Tx_codeword),cnt)
    end
end
save(resultsFileName,'params','BER','effectiveRate','Ensemble');

function [] = resultPlot(capacity,BER,numOfLines,n_block,rLDPC,rCC,overallR,cnt)
% plot group of deletions
lineType = ['-kd';'-go';'-b*';'-rs';'--c';'-xm'];
figure;
for idx = 1:numOfLines
    semilogy(capacity.',BER(:,numOfLines*((cnt - 1)/numOfLines - 1) + idx),lineType(mod(idx,length(lineType)),:),'LineWidth',2)
    hold on;
end
semilogy(overallR*ones(1,100),linspace(1e-4,1,100),'--k','LineWidth',1.5)
grid on; grid minor;
xlabel('C [nats]'); ylabel('BER')
legend('AWGN','BSC','BEC')
title(strcat('N_{block} = ',num2str(n_block),' R_{LDPC} = ',num2str(rLDPC),' R_{CC} = ',num2str(rCC),' R_{tot} = ',num2str(overallR)))
end