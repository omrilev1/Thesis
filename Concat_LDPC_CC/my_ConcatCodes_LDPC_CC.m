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

resultsFileName = strcat('.\results\',date,'TurboIterations_Results');


% Init Params
convCode.constraintLength = [8,7];
convCode.rate = [1/2];% [1/4 1/3 2/5 1/2];

decoding_method = [1]; % 0 for LDPC in hard decision , 1 for soft
dvb_rate = [2/5];%[9/10 5/6 3/4 2/3 2/5 ];

percentDeletions = [0 0.1 0.2 0.3];
turboIterations = [0,1,2,3];

blockIntrlvDepth = 256;

params = combvec(percentDeletions,convCode.constraintLength,convCode.rate,dvb_rate,decoding_method);
cnt = 1;

EbN0 = -0.5:0.25:4.5;
m = 1;

Nout = 50;
Nin = 100;

BER = zeros(length(EbN0),size(params,2));
BER_conv = zeros(length(EbN0),size(params,2));
effectiveRate = zeros(1,size(params,2));

NbitsLDPC = 64800;
while cnt<size(params,2)
    curr_del = params(1,cnt);
    curr_constraintLength = params(2,cnt);
    curr_conv_r = params(3,cnt);
    curr_LDPC_r = params(4,cnt);
    curr_decMethod = params(5,cnt);
    
    %% Generate the corresponding convolutional code
    [poly,tblLen] = convCodeGen(curr_constraintLength,curr_conv_r);
    
    % generate LDPC H matrix
    H = dvbs2ldpc(curr_LDPC_r);
    ldpc_enc = comm.LDPCEncoder(H);
    if curr_decMethod == 0
        ldpc_dec = comm.LDPCDecoder(H,'DecisionMethod','Hard decision','MaximumIterationCount',70);
    else
        ldpc_dec = comm.LDPCDecoder(H,'DecisionMethod','Soft decision','MaximumIterationCount',70); 
    end
    
    % generate convolutional code
    trellis = poly2trellis(curr_constraintLength+1,poly);
    conv_enc = comm.ConvolutionalEncoder('TrellisStructure',trellis);
    conv_dec = comm.APPDecoder(...
        'TrellisStructure',trellis, ...
        'Algorithm','Max*','CodedBitLLROutputPort',false);
    Nbits = size(H,2) - size(H,1);numOfErr_Conv = 0; numOfErr = 0;
    for i=1:length(EbN0)
        numOfErr_Conv = 0;
        numOfErr_Final = 0;
        for n_in = 1:Nin
            
            raw_bits = rand(Nbits,1) > 0.5;
            
            %% Encode LDPC
            ldpc_codedBits = step(ldpc_enc,raw_bits);
            
            %% Add deletions and interleave bits
            numOfDeletions = NbitsLDPC*curr_del;
            effectiveConvLen = NbitsLDPC - numOfDeletions;
            
            % truncate to use block interleaver - we add little more
            % deletions
            NbitsConv = blockIntrlvDepth*floor(effectiveConvLen/blockIntrlvDepth);
            numOfDeletions =  numOfDeletions + effectiveConvLen - NbitsConv ;
            
            % remove bits, that will be considered as deleted in the LDPC
            % decoding stage
            deletIdx = randperm(NbitsLDPC,numOfDeletions);
            dataIdx = setdiff(1:NbitsLDPC,deletIdx);
            
            ldpc_codedBits_ForConv = ldpc_codedBits(dataIdx);
            ldpc_codedBits_ForConv = (ldpc_codedBits_ForConv > 0.5);
            
            blockIntrlvMat = reshape(ldpc_codedBits_ForConv,[],blockIntrlvDepth);
            blockIntrlv    = blockIntrlvMat.';
            interleavedBits = blockIntrlv(:);
            
            %% Encode convolutional
%             Tx_codeword = convenc(interleavedBits,trellis);
            Tx_codeword = step(conv_enc,interleavedBits);
            
            %% channel
            x = qam_modulate(Tx_codeword(:),m,1);
            
            if m==1
                x = x.*exp(-1j*pi/4);
            end
            
            currR = Nbits/length(Tx_codeword);
            SNR = EbN0(i) + 10*log10(currR) + 10*log10(m);
            y = x + (1/sqrt(2))*10^(-SNR/20)*(randn(size(x)) + 1j*randn(size(x)));
            
            llr = calcLLR(y,m,10^(-SNR/10));
            
            %% Decode Convolutional code
            dataSoft = step(conv_dec,zeros(size(interleavedBits)),llr);
%             dataSoft = vitdec(llr,trellis,tblLen,'cont','unquant');
            dataHard = dataSoft > 0;
            if curr_decMethod == 0
                llr_ForDeintlrv = dataHard;
            else
                llr_ForDeintlrv = dataSoft;     
            end
            
            % count Convolutional code BER
            numOfErr_Conv = numOfErr_Conv + sum(dataHard(:) ~= interleavedBits(:));
            
            %% Deinterleave soft bits
            block_DeIntrlv = reshape(llr_ForDeintlrv,blockIntrlvDepth,[]);
            block_DeIntrlvMat = block_DeIntrlv.';
            LDPC_llrs = block_DeIntrlvMat(:);
            
            % add deletions - deletion will be treated as LLR = 0
            LDPC_final_llrs = zeros(NbitsLDPC,1);
            LDPC_final_llrs(dataIdx) = LDPC_llrs;
            LDPC_final_llrs(deletIdx) = 0;
            % decode LDPC
            ldpc_decodedLLRs = step(ldpc_dec,LDPC_final_llrs);
            
            ldpc_decodedBits = ldpc_decodedLLRs > 0;
            
            % count BER
            if curr_decMethod == 0
                numOfErr_Final = numOfErr_Final + sum(ldpc_decodedBits(:) ~= not(raw_bits));
            else
                numOfErr_Final = numOfErr_Final + sum(ldpc_decodedBits(:) ~= raw_bits);
            end
            if numOfErr_Final > 50
                BER(i,cnt) = numOfErr_Final/(Nbits*n_in);
                break
            end
        end
        if BER(i,cnt) < 1e-6
            BER(i:end,cnt) = 1e-6;
            BER_conv(i,cnt) = numOfErr_Conv/(n_in*length(interleavedBits));
            
            effectiveRate(cnt) = Nbits/length(Tx_codeword);
            iterMessage = sprintf('LDPC : rate %0.2f , Deletion Percentage = %0.2f \nConvolutional Rate = %0.2f constraint length = %0.2f \n Effective Rate = %0.2f%',curr_LDPC_r,100*curr_del,curr_conv_r,curr_constraintLength,effectiveRate(cnt));
            iterMessage2 = sprintf('BER Convolutional = %.5f overall BER = %.5f',BER_conv(i,cnt),BER(i,cnt));
            iterMessage3 = sprintf('EbN0 = %2g',EbN0(i));
            disp(iterMessage3);disp(iterMessage);disp(iterMessage2);disp('%%%%%%%%%%%%%%%%%%%%%');
            break
            
        end
        
        BER_conv(i,cnt) = numOfErr_Conv/(n_in*length(interleavedBits));
        effectiveRate(cnt) = Nbits/length(Tx_codeword);
        iterMessage = sprintf('LDPC : rate %0.2f , Deletion Percentage = %0.2f \nConvolutional Rate = %0.2f constraint length = %0.2f \n Effective Rate = %0.2f%',curr_LDPC_r,100*curr_del,curr_conv_r,curr_constraintLength,effectiveRate(cnt));
        iterMessage2 = sprintf('BER Convolutional = %.5f overall BER = %.5f',BER_conv(i,cnt),BER(i,cnt));
        iterMessage3 = sprintf('EbN0 = %2g',EbN0(i));
        disp(iterMessage3);disp(iterMessage);disp(iterMessage2);disp('%%%%%%%%%%%%%%%%%%%%%');
    end
    cnt = cnt + 1;
    
    if mod(cnt-1,4) == 0
    % plot group of deletions
    lineType = ['-kd';'-go';'-b*';'-rs';'--c'];
    figure;
    for idx = 1:(length(lineType) - 1)
        semilogy(EbN0.',BER(:,4*((cnt - 1)/4 - 1) + idx),lineType(idx,:),'LineWidth',2)
        hold on;
    end
    semilogy(EbN0.',BER_conv(:,(cnt - 1)),'-xm','LineWidth',2)
    grid on; grid minor;
    xlabel('E_b / N_0 [dB]'); ylabel('BER')
    
    legend(sprintf('Deletion Percentage = %0.2f  , effective Rate = %0.2f',percentDeletions(1),effectiveRate(4*((cnt - 1)/4 - 1) + 1)),...
        sprintf('Deletion Percentage = %0.2f  , effective Rate = %0.2f',percentDeletions(2),effectiveRate(4*((cnt - 1)/4 - 1) + 2)),...
        sprintf('Deletion Percentage = %0.2f  , effective Rate = %0.2f',percentDeletions(3),effectiveRate(4*((cnt - 1)/4 - 1) + 3)),...
        sprintf('Deletion Percentage = %0.2f  , effective Rate = %0.2f',percentDeletions(4),effectiveRate(4*((cnt - 1)/4 - 1) + 4)),...
        'convolutional BER');
    
    plotMessage = sprintf('LDPC : rate %0.2f \nConvolutional Rate = %0.2f \n constraint length = %0.2f',curr_LDPC_r,curr_conv_r,curr_constraintLength);
    title(plotMessage)
    end
end
save(resultsFileName,'params','BER','effectiveRate');

