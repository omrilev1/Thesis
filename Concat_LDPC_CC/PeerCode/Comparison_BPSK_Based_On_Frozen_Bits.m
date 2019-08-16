close all;
%clear all;

%%
%----------------------------Parallel Computing----------------------------%

myCluster = parcluster('local');
NumOfWorkers = myCluster.NumWorkers;
try
    matlabpool local
end

%%
%--------------------------------N = 2^10 , 16QAM--------------------------------%

%%
%N = 2^10; % Number of bits in block
% Constellation_Type = 'QAM';
%m = 4;
%R_Channel_Bits = 4*(1/8); % [bits/channel use]
%SNR_Vec_dB_EbN0 = [1,2,3,4,5,6];
%
%BLER_Max = 10^(-5);

%%
% N = 2^10; % Number of bits in block
% Constellation_Type = 'QAM';
% m = 4;
% R_Channel_Bits = m*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [3,4,5];
% 
% BLER_Max = 10^(-5);

%%
% N = 2^10; % Number of bits in block
% Constellation_Type = 'QAM';
% m = 4;
% R_Channel_Bits = m*(7/8); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [6.5,7.5,8.5,9.5,10.5];
% BLER_Max = 10^(-5);

%%
%--------------------------------N = 2^8 , 16QAM--------------------------------%
% N = 2^8; % Number of bits in block
% Constellation_Type = 'QAM';
% m = 4;
% R_Channel_Bits = m*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [3,4,5];
% 
% BLER_Max = 10^(-5);

%%
%--------------------------------N = 2^10 , 256QAM------------------------------%

%%
%--------------------------------N = 2^8 , 256QAM-------------------------------%
% N = 2^8; % Number of bits in block
% Constellation_Type = 'QAM';
% m = 8;
% R_Channel_Bits = m*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [9,10,11];
% 
% BLER_Max = 10^(-5);

%%
%--------------------------------N = 2^10 , 16PAM--------------------------------%

%%
% N = 2^10; % Number of bits in block
% Constellation_Type = 'PAM';
% m = 4;
% R_Channel_Bits = m*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [9,10,11];
% 
% BLER_Max = 10^(-5);

%%
% N = 2^10; % Number of bits in block
% Constellation_Type = 'PAM';
% m = 4;
% R_Channel_Bits = m*(3/4); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [13,14,15];
% 
% BLER_Max = 10^(-5);

%%
%Try
N = 2^10; % Number of bits in block
Constellation_Type = 'QAM';
m = 4;
R_Channel_Bits = m*(1/2); % [bits/channel use]
SNR_Vec_dB_EbN0 = [4];

BLER_Max = 10^(-5);

m_BIPCM = m/2;

%%

R_Code_Bits = R_Channel_Bits/m; % [Information bits / Total bits]

SNR_Vec_dB_EsN0 = SNR_Vec_dB_EbN0 + 10*log10(m) + 10*log10(R_Code_Bits);

%%
SNR_Vec = 10.^(SNR_Vec_dB_EsN0./10);

BLER_Axis_Buffer = 10^-1;
Min_BLER = BLER_Axis_Buffer*min(BLER_Max);

Min_Num_Of_Errors_Sufficient_Statistics = 100;

file_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Files');
plot_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Plots');
code_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Code');
addpath(genpath(code_path));

is_visible = 'on';

BLER = inf(2,length(SNR_Vec));
BER = inf(2,length(SNR_Vec));

SNR_Type = 'Eb/N0';

date_string = strrep(sprintf('%s',datestr(now)),':','_');
file_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%d%s_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),2^m,Constellation_Type,date_string);
plot_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%d%s_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),2^m,Constellation_Type,date_string);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
%Try

[BLER(1,:),BER(1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,Constellation_Type,m,zeros(0),zeros(0),'SP',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 

load(fullfile(file_path,sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R_Code_Bits),'.','_'),2^m,Constellation_Type,'SP')),'Is_Frozen_Bit_Index_Vec');
Is_Frozen_Bit_Index_Vec_QAM = reshape(Is_Frozen_Bit_Index_Vec,N/m,[]).';
R_QAM = sum(~Is_Frozen_Bit_Index_Vec_QAM,2)/(N/m);

Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,'SP',zeros(0));
d_BPSK = 2;

for i=1:1:m
    d_QAM(i) = abs(Constellation_Mapping_Array(1+0)-Constellation_Mapping_Array(1+2^(m-i)));
    [BLER(1+i,:),BER(1+i,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'gArikan',1,1,0,'BPSK',1,zeros(0),zeros(0),'BPSK',N/m,R_QAM(i),SNR_Vec.*((d_QAM(i)/d_BPSK).^2),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
    save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 
end


%%
%----------------------------plot BLER----------------------------%

legend_name = {sprintf('MLPC %d%s SP (regular frozen bits)',2^m,Constellation_Type),sprintf('MLPC %d%s SP (BPSK based on frozen bits)',2^m,Constellation_Type)};
title_name = sprintf('N=%d[bits] R=%s[information bits per channel use]',N,sprintf('%.3f',R_Channel_Bits));

save(fullfile(file_path,file_name),'plot_name','legend_name','title_name','N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,length(BLER(:,1)),Min_BLER,title_name,legend_name,plot_path,plot_name);

%%
sprintf('finish')
    