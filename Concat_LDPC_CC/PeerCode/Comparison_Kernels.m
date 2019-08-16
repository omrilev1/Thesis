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
N = 2*4^5; % Number of bits in block
% N = 2*4^4; % Number of bits in block
% Decode_Type = 'BIMLPCM_gArikan_GF4';
% Lable_Type = 'SP-Compound-GF4';
% Decode_Type = 'MLPC_gArikan_GF4';
% Lable_Type = 'SP-GF4';
% Decode_Type = 'BIPCM_gArikan_GF4';
% Lable_Type = 'Gray-GF4';
% Constellation_Type = 'QAM';
% m = 8;
% Constellation_Type = 'PAM';
% m = 8;
% Constellation_Type = 'QAM';
% m = 4;
% Constellation_Type = 'PAM';
% m = 4;
Constellation_Type = 'BPSK';
m = 1;
R_Channel_Bits = m*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [6,6.5,7,7.5];
% SNR_Vec_dB_EbN0 = [9,10,11];
SNR_Vec_dB_EbN0 = [1.75,2];
% SNR_Vec_dB_EbN0 = [8,8.5,9,9.5,10,10.5,11,11.5,12,12.5];
% SNR_Vec_dB_EbN0 = [3.5];
% SNR_Vec_dB_EbN0 = [8,8.5,9,9.5,10,10.5];

BLER_Max = 10^(-5);

m_BIPCM = m/2;
% m_BIPCM = 2;

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

Num_Of_Labels_To_Compare = 4;

BLER = inf(Num_Of_Labels_To_Compare,length(SNR_Vec));
BER = inf(Num_Of_Labels_To_Compare,length(SNR_Vec));

SNR_Type = 'Eb/N0';


date_string = strrep(sprintf('%s',datestr(now)),':','_');
file_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_%d%s_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),2^m,Constellation_Type,date_string);
plot_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_%d%s_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),2^m,Constellation_Type,date_string);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
%Try

% [BLER(1,:),BER(1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'gArikan',1,1,0,Constellation_Type,m,zeros(0),zeros(0),'Gray',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 
% 
% [BLER(2,:),BER(2,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'gRS4',1,1,0,Constellation_Type,m,zeros(0),zeros(0),'Gray',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 
% 
% [BLER(3,:),BER(3,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'gArikan_GF4',1,1,0,Constellation_Type,m,zeros(0),zeros(0),'Gray',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 

[BLER(1,1:3),BER(1,4:6)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'BIPCM',1,1,0,Constellation_Type,m,zeros(0),zeros(0),'Gray',N,R_Code_Bits,SNR_Vec(4:6),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 

[BLER(2,1:3),BER(2,4:6)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'BIPCM_gArikan_GF4',1,1,0,Constellation_Type,m,zeros(0),zeros(0),'Gray',N,R_Code_Bits,SNR_Vec(4:6),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 

% [BLER(1,:),BER(1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'BIPCM_gRS4',1,1,0,Constellation_Type,m,zeros(0),zeros(0),'Gray',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 

[BLER(3,4:6),BER(3,1:3)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,Constellation_Type,m,zeros(0),zeros(0),'SP',N,R_Code_Bits,SNR_Vec(1:3),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 

[BLER(4,4:6),BER(4,1:3)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC_gArikan_GF4',1,1,0,Constellation_Type,m,zeros(0),zeros(0),'SP',N,R_Code_Bits,SNR_Vec(1:3),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 

% [BLER(2,:),BER(2,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC_gRS4',1,1,0,Constellation_Type,m,zeros(0),zeros(0),'SP',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 

% [BLER(1,:),BER(1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'BIMLPCM',1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),'SP-Compound',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 
% 
% [BLER(3,:),BER(3,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'BIMLPCM_gRS4',1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),'SP-Compound',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 
% 
% [BLER(3,:),BER(3,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'BIMLPCM_gArikan_GF4',1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),'SP-Compound',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 

% [BLER(3,:),BER(3,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'BIMLPCM_gRS4',1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),'SP-Gray-GF4',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 


%%
%----------------------------plot BLER----------------------------%

% legend_name = {sprintf('gArikan'),sprintf('gRS4'),sprintf('gArikan GF4')};
% legend_name = {sprintf('BIPCM UVV %d%s Gray', 2^m,Constellation_Type),sprintf('BIPCM UVV4 %d%s Gray', 2^m,Constellation_Type),sprintf('BIPCM RS4 %d%s Gray', 2^m,Constellation_Type),sprintf('MLPC UVV %d%s SP', 2^m,Constellation_Type),sprintf('MLPC UVV4 %d%s SP', 2^m,Constellation_Type),sprintf('MLPC RS4 %d%s SP', 2^m,Constellation_Type)};
legend_name = {sprintf('BIPCM UVV %d%s Gray', 2^m,Constellation_Type),sprintf('BIPCM UVV4 %d%s Gray', 2^m,Constellation_Type),sprintf('MLPC UVV %d%s SP', 2^m,Constellation_Type),sprintf('MLPC UVV4 %d%s SP', 2^m,Constellation_Type)};
title_name = sprintf('N=%d[bits] R=%s[information bits per channel use]',N,sprintf('%.3f',R_Channel_Bits));

save(fullfile(file_path,file_name),'plot_name','legend_name','title_name','N','R_Code_Bits','Constellation_Type','m','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

% Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,length(BLER(:,1)),Min_BLER,title_name,legend_name,plot_path,plot_name);
Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,2,Min_BLER,title_name,legend_name,plot_path,plot_name);

%%
sprintf('finish')
    