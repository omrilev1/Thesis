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
%Try
N = 2^10; % Number of bits in block

m = 4;

R_Channel_Bits = m*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [1,1.5,2,2.5,3,3.5,4,4.5,5];
% SNR_Vec_dB_EbN0 = [2,2.5,3.5,4,4.5,5,5.5];
% SNR_Vec_dB_EbN0 = [-1.5,-1,-0.5,0];
SNR_Vec_dB_EbN0 = [1.5];
%%

R_Code_Bits = R_Channel_Bits/m; % [Information bits / Total bits]

SNR_Vec_dB_EsN0 = SNR_Vec_dB_EbN0 + 10*log10(m) + 10*log10(R_Code_Bits);

%%
SNR_Vec = 10.^(SNR_Vec_dB_EsN0./10);

BLER_Axis_Buffer = 10^-1;
BLER_Max = 10^(-5);
Min_BLER = BLER_Axis_Buffer*min(BLER_Max);

Min_Num_Of_Errors_Sufficient_Statistics = 100;

file_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Files');
plot_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Plots');
code_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Code');
addpath(genpath(code_path));

% is_visible = 'off';
is_visible = 'on';

m_BIPCM = m/2;

Num_Of_Labels_To_Compare = 3;

BLER = inf(Num_Of_Labels_To_Compare,length(SNR_Vec));
BER = inf(Num_Of_Labels_To_Compare,length(SNR_Vec));

SNR_Type = 'Eb/N0';

date_string = strrep(sprintf('%s',datestr(now)),':','_');
file_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%d-bit_Multi_Dimensional_Constellations_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),m,date_string);
plot_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%d-bit_Multi_Dimensional_Constellations_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),m,date_string);
save(fullfile(file_path,file_name),'N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
%Try

[BLER(1,:),BER(1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'Multi_Dimensional_BIMLPCM',1,1,0,sprintf('%d-bit-2x8PSK',m),m,m_BIPCM,zeros(0),'Multi-Dimensional-SP-Compound',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');



[BLER(1,:),BER(1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'BIPCM',1,1,0,'QAM',2,m_BIPCM,zeros(0),'Gray',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

[BLER(2,:),BER(2,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,'QAM',2,m_BIPCM,zeros(0),'SP',N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%

% [BLER(1,6:9),BER(1,6:9)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'BIPCM',1,1,0,'QAM',m,m_BIPCM,zeros(0),'Gray',N,R_Code_Bits,SNR_Vec(6:9),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');
% 
% [BLER(2,2:5),BER(2,2:5)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'Multi_Dimensional_BIPCM',1,1,0,sprintf('%d-bit-2x8PSK',m),m,m_BIPCM,zeros(0),'Multi-Dimensional-Gray',N,R_Code_Bits,SNR_Vec(2:5),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');
% 
% [BLER(3,1:4),BER(3,1:4)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'Multi_Dimensional_BIPCM',1,1,0,sprintf('%d-bit-3x4PSK',m),m,m_BIPCM,zeros(0),'Multi-Dimensional-Gray',N,R_Code_Bits,SNR_Vec(1:4),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

[BLER(1,3:6),BER(1,3:6)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'BIPCM',1,1,0,'QAM',m,m_BIPCM,zeros(0),'Gray',N,R_Code_Bits,SNR_Vec(3:6),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

[BLER(2,3:6),BER(2,3:6)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'Multi_Dimensional_BIPCM',1,1,0,sprintf('%d-bit-2x8PSK',m),m,m_BIPCM,zeros(0),'Multi-Dimensional-Gray',N,R_Code_Bits,SNR_Vec(3:6),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

[BLER(3,1:4),BER(3,1:4)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'Multi_Dimensional_BIPCM',1,1,0,sprintf('%d-bit-3x4PSK',m),m,m_BIPCM,zeros(0),'Multi-Dimensional-Gray',N,R_Code_Bits,SNR_Vec(1:4),BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
%----------------------------plot BLER----------------------------%

legend_name{Num_Of_Labels_To_Compare} = [];

legend_name{1} = sprintf('BIPCM %dQAM Gray',2^m);
legend_name{2} = sprintf('BIPCM %d-bit-2x8PSK Gray',m);
legend_name{3} = sprintf('BIPCM %d-bit-3xQPSK Gray',m);

title_name = sprintf('N=%d[bits] R=%s[information bits per channel use]',N,sprintf('%.3f',R_Channel_Bits));

save(fullfile(file_path,file_name),'plot_name','legend_name','title_name','N','R_Code_Bits','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER'); 

Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,length(BLER(:,1)),Min_BLER,title_name,legend_name,plot_path,plot_name);

%%
sprintf('finish')
    