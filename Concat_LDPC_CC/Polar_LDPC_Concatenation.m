clear all; clc;
addpath(genpath('.\QAM'))
addpath(genpath('.\results'))
addpath(genpath('.\Ensembles'))
addpath(genpath('.\PeerCode'))

%% This script simulate concatenated coding scheme :
% The outer code is high rate LDPC. The inner code suppose to be
% some code with performance guarantee over BSC (for example convolutional codes or polar codes).
% The main thing we check is the universality of the code

% The overall coding scheme is linear, thus we simulate only the zeros codeword
% and randomize the noise

% The simulation can be done for AWGN, BSC and BEC

resultsFileName = strcat('.\results\Polar',date,'GallagerEnsemble');

%%
%----------------------------Parallel Computing----------------------------%
toUseParallel = 1;

if toUseParallel
    myCluster = parcluster('local');
    NumOfWorkers = myCluster.NumWorkers;
    try
        parpool local
    end
else
    NumOfWorkers = 1;
end

Rtot = 0.5;

%%  Init Params
%% LDPC parameters
fullBP_iterations = 0; % number of iterations between inner and outer codes
decoding_method = 1; % 0 for LDPC in hard decision , 1 for soft

% LDPC Ensemble
Ensemble = 'Irregular_BEC';% 'Irregular_BEC';% 'Gallager';
ldpc_rate = 0.66;% [0.87,0.9,0.82,0.66];

puncture_percent = 0.19;% [0.22 0.15 0.1];

%% Polar Code parameters
% list decoder related parameters
L_List = 32;
CRC_length = 0;
L_Optimize_List = ones(1,size(L_List,2));
GA_CRC_Length_List = CRC_length*ones(1,size(L_List,2)); %crc16

% code related parameters
N_polar = 2^11; % Number of bits in block
Constellation_Type = 'BPSK';
Decode_Type = 'gArikan';
Lable_Type = 'BPSK';
m = 1;
R = 0.55; % desired code rate
RtotPolar = (R*N_polar - CRC_length)/N_polar; % total code rate after CRC addition
R_Channel_Bits = m*R; % [bits/channel use]
m_SCL = 1;
R_Code_Bits = R_Channel_Bits/m; % [Information bits / Total bits]

U_N = N_polar;
Is_GF4 = false;
g_Arikan = [1,0;1,1];
M = 2^m;

% Generate polar encoding matrix
g = g_Arikan;
for i=1:1:log2(N_polar)-1
    g = kron(g,g_Arikan);
end
g = bitrevorder(g);

g0 = 1;
Constellation_Mapping_Array = [-1,1];

%% Concatenation parameters
delta1 = 0.5736; % lower Bahattacharya threshold
delta2 = 0.83; % higher Bahattacharya threshold ( delta1 < Z < delta2 ==> protect with LDPC)

%% channel parameters
channelType = [2,1,0]; % {0 == 'BEC',1 == 'BSC',2 == 'AWGN'};
cap = 0.46:0.0125:0.75;
blockIntrlvDepth = 128;
params = combvec(fullBP_iterations,puncture_percent,decoding_method,channelType,ldpc_rate);

cnt = 1;

m = 1;

Nout = 5;
Nin = 500;

BLER_Max = 10^(-2);

BLER_Axis_Buffer = 10^-1;
Min_BLER = BLER_Axis_Buffer*min(BLER_Max);

Min_Num_Of_Errors_Sufficient_Statistics = 100;

BER_total = zeros(length(cap),size(params,2));
BER_polar = zeros(length(cap),size(params,2));
effectiveRate = zeros(1,size(params,2));

while cnt<=size(params,2)
    curr_iter = params(1,cnt);
    curr_punct = params(2,cnt);
    curr_decMethod = params(3,cnt);
    curr_ldpc_rate = params(5,cnt);
    curr_chan = params(4,cnt);
    
    curr_chan_string = createChannelString(curr_chan);
    %% Generate the corresponding convolutional code
    % generate LDPC H matrix
    [H,curr_LDPC_r,ldpc_dec] = ldpcGen(curr_ldpc_rate,Ensemble,curr_decMethod);
    
    N_LDPC = size(H,2);numOfErr_Polar = 0; numOfErr = 0;
    Ndata = N_LDPC - size(H,1);
    
    SER_Vec = zeros(size(cap,2),U_N);
    Is_Frozen_Bit_Index_Vec = zeros(size(cap,2),N_polar);
    Is_Frozen_Bit_Index_Vec = Is_Frozen_Bit_Index_Vec == ones(size(cap,2),N_polar);
    Is_Frozen_Bit_Index_List_Vec = zeros(size(L_List,2),size(cap,2),N_polar);
    Is_Frozen_Bit_Index_List_Diff_Vec = zeros(size(L_List,2),size(cap,2));
    Is_Frozen_Bit_Index_temp_Vec = zeros(1,N_polar);
    Genie_Aided_LLR_Avg_Vec = zeros(size(cap,2),U_N);
    Genie_Aided_LLR_Var_Vec = zeros(size(cap,2),U_N);
    Estimated_SER_Vec = zeros(size(cap,2),U_N);
    BLER = zeros(size(L_List,2),size(cap,2));
    BER = zeros(size(L_List,2),size(cap,2));
    First_Symbol = zeros(size(L_List,2),size(cap,2),U_N);
    temp_First_Symbol = zeros(1,U_N);
    temp_LLR_Avg = zeros(1,U_N);
    temp_LLR_Var = zeros(1,U_N);
    LLR_Avg = zeros(size(L_List,2),size(cap,2),U_N);
    LLR_Var = zeros(size(L_List,2),size(cap,2),U_N);
    Is_Change_Frozen_Bits_Index_Vec = zeros(size(cap,2)-1,N_polar);
    Estimated_Bhattacharyya = zeros(size(cap,2),U_N);
    is_high_SNR = zeros(1,size(cap,2));
    
    for i=1:length(cap)
        numOfErr_Polar = zeros(1,curr_iter+1);
        numOfErr_Final = zeros(1,curr_iter+1);
        
        %% Find optimal set of frozen bits for current run
        [Estimated_Bhattacharyya(i,:),SER_Vec(i,:),Is_Frozen_Bit_Index_Vec(i,:),~,Genie_Aided_LLR_Avg_Vec(i,:),Genie_Aided_LLR_Var_Vec(i,:),...
            Estimated_SER_Vec(i,:),~,~,~,~,~] = listDecoderSim_GeneralChannel(NumOfWorkers,true,g,g0,Decode_Type,1,0,Constellation_Type,m,[],m_SCL,Constellation_Mapping_Array,N_polar,R,BLER_Max,...
            cap(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_Vec(i,:),curr_chan);
        
        Is_Frozen_Bit_Index_Vec(i,:) = Is_Frozen_Bit_Index_Vec(i,:) == ones(1,N_polar);
        %     LDPC_idx = find((Estimated_Bhattacharyya(i,:) > delta1) & (Estimated_Bhattacharyya(i,:) < delta2) ...
        %                     & (~Is_Frozen_Bit_Index_Vec(i,:)));
        %     LDPC_idx = (~Is_Frozen_Bit_Index_Vec(i,:));
        %     polar_Idx = (~LDPC_idx) & (~Is_Frozen_Bit_Index_Vec(i,:));
        
        %% Rate Matching (Puncture LDPC and take only the relevant bits into the polar)
        % (Note that we can ue lower LDPC code with longr block
        % length)
        numOfPunct = N_LDPC - sum((~Is_Frozen_Bit_Index_Vec(i,:)));
        
        % remove bits, that will be considered as punctured in the LDPC
        % decoding stage
        punctIdx = randperm(N_LDPC,numOfPunct);
        dataIdx = setdiff(1:N_LDPC,punctIdx);
        
        k_polar = sum((~Is_Frozen_Bit_Index_Vec(i,:)));
        ilvIdx = randperm(k_polar);
        deIlv(ilvIdx) = 1:k_polar;
        for n_out = 1:Nout
            for n_in = 1:Nin
                
                %% Encode LDPC
                ldpc_codedBits = zeros(1,N_LDPC);
                
                ldpc_codedBits_ForPolar = ldpc_codedBits(dataIdx);
                ldpc_codedBits_ForPolar = (ldpc_codedBits_ForPolar > 0.5);
                
                interleavedBits = ldpc_codedBits_ForPolar(ilvIdx);
                
                %% Encode Polar
                % simple scrambler
                scramblingSeq = rand(size(interleavedBits)) > 0.5;
                scrambledBits = mod(scramblingSeq + interleavedBits,2);
                
                ui = rand(1,N_polar) > 0.5;
                ui(~Is_Frozen_Bit_Index_Vec(i,:)) = scrambledBits;
                
                Tx_codeword = mod(ui*g,2);
                
                %% Modulate and send through the channel
                
                x = (-1).^(Tx_codeword(:));
                
                %% channel
                
                currR = Ndata/length(Tx_codeword);
                [llr,ll0,ll1] = channelModel(x,cap(i),curr_chan);
                
                %% Decode Polar + LDPC in an iterative manner
                
                polarCode_priors = zeros(size(interleavedBits));
                for decodeIter=0:curr_iter
                    
                    % Decode Polar code
                    
                    % in the first iteraticon take aprior probabilities to be
                    % [1/2,1/2]. After that take the apriori probabilities from
                    % LDPC output LLRs
                    
                    [dataHard,~,dataSoft] = gArikan_SC_Decoder(llr,ui,N_polar,Is_Frozen_Bit_Index_Vec(i,:),0);
                    
                    % take only information part
                    if curr_decMethod == 0
                        llr_ForDeintlrv = dataHard(~Is_Frozen_Bit_Index_Vec(i,:));
                    else
                        llr_ForDeintlrv = dataSoft(~Is_Frozen_Bit_Index_Vec(i,:));
                    end
                    llr_ForDeintlrv = llr_ForDeintlrv.*(2*scramblingSeq - 1);
                    
                    % count Convolutional code BER
                    numOfErr_Polar(decodeIter+1) = numOfErr_Polar(decodeIter+1) + ...
                        sum(dataHard(~Is_Frozen_Bit_Index_Vec(i,:)) ~= scrambledBits);
                    
                    %% Deinterleave soft bits
                    LDPC_llrs = -1*llr_ForDeintlrv(deIlv);
                    
                    %% Rx punctures - treated as LLR = 0
                    LDPC_final_llrs = zeros(N_LDPC,1);
                    LDPC_final_llrs(dataIdx) = LDPC_llrs;
                    LDPC_final_llrs(punctIdx) = -0.001;
                    
                    % decode LDPC
                    [ldpc_decodedLLRs,numIterations] = step(ldpc_dec,LDPC_final_llrs);
                    ldpc_decodedBits = ldpc_decodedLLRs < 0;
                    
                    % count BER
                    if (sum(abs(ldpc_decodedLLRs - LDPC_final_llrs)) == 0)
                        numOfErr_Final(decodeIter+1) = numOfErr_Final(decodeIter+1) + Ndata/2;
                    else
                        numOfErr_Final(decodeIter+1) = numOfErr_Final(decodeIter+1) + sum(ldpc_decodedBits(2:Ndata-1));
                    end
                    
                    % calculate polar code next iteration llrs
                    polarCode_priors = ldpc_decodedLLRs(dataIdx);
                    polarCode_priors = polarCode_priors(ilvIdx);
                end
            end
            if numOfErr_Final(end) > 50
                BER_total(i,cnt) = numOfErr_Final(end)/((Ndata-2)*Nin);
                break
            end
        end
        
        BER_polar(i,cnt) = numOfErr_Polar(decodeIter+1)/(n_in*n_out*length(interleavedBits));
        effectiveRate(cnt) = Ndata/length(Tx_codeword);
        iterMessage1 = sprintf('LDPC : rate %0.2f , Punct = %0.2f \n Polar : Rate = %0.2f N_{polar} = %0.2f \n Effective Rate = %0.2f%',...
            curr_LDPC_r,100*curr_punct,R,N_polar,effectiveRate(cnt));
        iterMessage2 = sprintf('BER Polar = %.5f overall BER = %.5f',BER_polar(i,cnt),BER_total(i,cnt));
        iterMessage3 = sprintf('C = %2g',cap(i));
        iterMessage4 = sprintf(strcat('Channel is ',createChannelString(curr_chan)));
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

function [chan_string] = createChannelString(chan)

if chan == 0
    chan_string = 'BEC';
elseif chan == 1
    chan_string = 'BSC';
else
    chan_string = 'AWGN';
end
end