close all; clear all; clc;
addpath(genpath('.\QAM'))
addpath(genpath('.\results'))


% Init Params
constraintLength = [5,6,8,10,12,14];
rsc = [0];
rate =  [1/2];% [1/4 1/3 2/5 1/2];

params = combvec(constraintLength,rate,rsc);
cnt = 1;

% EbN0 = (2:0.5:5.5);
m = 1;
cap = 0.45:0.025:0.65;
Nout = 50;
Nin = 5;
curr_chan = 1;
% BER = zeros(length(EbN0),size(params,2));
BER = zeros(length(cap),size(params,2));

downsamplePlot = 6;
Nbits = 2^14; % 64800
while cnt<=size(params,2)
    curr_constraintLength = params(1,cnt);
    curr_conv_r = params(2,cnt);
    curr_conv_rsc = params(3,cnt);
    %% Generate the corresponding convolutional code
    [trellis,conv_enc,conv_dec ] = convCodeGen(curr_constraintLength,curr_conv_r,curr_conv_rsc);
    
%     for i=1:length(EbN0)
    for i=1:length(cap)
        numOfErr_Conv = 0;
        for n_out = 1:Nout
            for n_in = 1:Nin
                
                raw_bits = randi([0 1],Nbits,1);
                
                %% Encode convolutional
                Tx_codeword = step(conv_enc,raw_bits);
                
                %% channel
                x = qammod(Tx_codeword.',2,'InputType','bit','UnitAveragePower',true);
                
                currR = Nbits/length(Tx_codeword);
                
                [llr] = channelModel(x,cap(i),curr_chan);
                
%                 SNR = EbN0(i) + 10*log10(currR) + 10*log10(m);
%                 y = awgn(x,SNR,'measured');
%                 llr = qamdemod(-1*y,2,'OutputType','approxllr', ...
%                     'UnitAveragePower',true,'NoiseVariance',10^(-SNR/10));
                %% Decode Convolutional code + LDPC in an iterative manner
                
                convCode_priors = zeros(size(raw_bits));
                
%                 dataSoft = step(conv_dec,convCode_priors,llr.');
                dataSoft = step(conv_dec,convCode_priors,llr);
                
                dataHard = dataSoft > 0;
                
                % count BER
                numOfErr_Conv = numOfErr_Conv + sum(dataHard(10:end) ~= raw_bits(10:end));
            end
            if numOfErr_Conv > 20
                BER(i,cnt) = numOfErr_Conv/((Nbits-10)*n_in*n_out);
                break
            end
        end
        
        BER(i,cnt) = numOfErr_Conv/(n_in*n_out*(length(raw_bits)-10));
%         iterMessage = sprintf('rsc %d rate %0.2f , v = %0.2f , EbN0 = %0.2f [dB] BER = %0.6f',...
%             curr_conv_rsc,curr_conv_r,curr_constraintLength,EbN0(i),BER(i,cnt));
        iterMessage = sprintf('rsc %d rate %0.2f , v = %0.2f , C = %0.2f [bits/channel use] BER = %0.6f',...
            curr_conv_rsc,curr_conv_r,curr_constraintLength,cap(i),BER(i,cnt));
        disp(iterMessage);disp('%%%%%%%%%%%%%%%%%%%%%');
    end
    if mod(cnt,downsamplePlot) == 0
        % plot group of deletions
        lineType = ['-kd';'-go';'-b*';'--c';'-.r';'-ms'];
        figure;
%         for idx = 1:downsamplePlot
%             semilogy(EbN0.',BER(:,downsamplePlot*(cnt/downsamplePlot - 1) + idx),lineType(idx,:),'LineWidth',2)
%             hold on;
%         end
%         grid on; grid minor;
%         xlabel('E_b / N_0 [dB]'); ylabel('BER')
        for idx = 1:downsamplePlot
            semilogy(cap.',BER(:,downsamplePlot*(cnt/downsamplePlot - 1) + idx),lineType(idx,:),'LineWidth',2)
            hold on;
        end
        grid on; grid minor;
        xlabel('C [bits/channel use]'); ylabel('BER')        
        legend('v = 5','v = 6','v = 8','v = 10','v = 12','v = 14');
        title(strcat('rate= ',num2str(curr_conv_r),' n = ',num2str(Nbits)))
        
        figure;
%         for idx = 1:downsamplePlot
%             semilogy(EbN0.',2*sqrt(BER(:,downsamplePlot*(cnt/downsamplePlot - 1) + idx).*(1-BER(:,downsamplePlot*(cnt/downsamplePlot - 1) + idx))),...
%                 lineType(idx,:),'LineWidth',2)
%             hold on;
%         end
%         grid on; grid minor;
%         xlabel('E_b / N_0 [dB]'); ylabel('Equivalent Bhattacharya')
        for idx = 1:downsamplePlot
            semilogy(cap.',2*sqrt(BER(:,downsamplePlot*(cnt/downsamplePlot - 1) + idx).*(1-BER(:,downsamplePlot*(cnt/downsamplePlot - 1) + idx))),...
                lineType(idx,:),'LineWidth',2)
            hold on;
        end
        grid on; grid minor;
        xlabel('C [bits/channel use]'); ylabel('Equivalent Bhattacharya')        
        legend('v = 5','v = 6','v = 8','v = 10','v = 12','v = 14');
        title(strcat('rate= ',num2str(curr_conv_r),' n = ',num2str(Nbits)))
    end
    cnt = cnt + 1;
end


