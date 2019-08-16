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
%--------------------------------N = 2^11--------------------------------%
% N = 2^11; % Number of bits in block
% m_SNR_EbN0 = 1;
% R_Channel_Bits = m_SNR_EbN0*(1/4); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [1.25,2,2.75];
% m_SNR_EbN0 = 1;

%%
%--------------------------------N = 2^11--------------------------------%
% N = 2^11; % Number of bits in block
% m_SNR_EbN0 = 1;
% R_Channel_Bits = m_SNR_EbN0*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [1.25,2,2.75];
% m_SNR_EbN0 = 1;

%%
%--------------------------------N = 2^11--------------------------------%
% N = 2^11; % Number of bits in block
% m_SNR_EbN0 = 1;
% R_Channel_Bits = m_SNR_EbN0*(3/4); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [2.5,3,3.5];

%%
%--------------------------------N = 2^11--------------------------------%
% N = 2^11; % Number of bits in block
% m_SNR_EbN0 = 2;
% R_Channel_Bits = m_SNR_EbN0*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [-1.5,-0.75,0];

%%
%--------------------------------N = 2^11--------------------------------%
% N = 2^11; % Number of bits in block
% m_SNR_EbN0 = 2;
% R_Channel_Bits = m_SNR_EbN0*(3/4); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [-0.5,0.5,1.5];

%%
%Try
N = 2^11; % Number of bits in block
m_SNR_EbN0 = 1;
R_Channel_Bits = m_SNR_EbN0*(1/2); % [bits/channel use]
SNR_Vec_dB_EbN0 = [1.75,2,2.25];
m_SNR_EbN0 = 1;
%%

SNR_Vec_dB_EsN0 = SNR_Vec_dB_EbN0 + 10*log10(m_SNR_EbN0) + 10*log10(R_Channel_Bits);
SNR_Vec = 10.^(SNR_Vec_dB_EsN0./10);
BLER_Max = 10^(-5);

BLER_Axis_Buffer = 10^-1;
Min_BLER = BLER_Axis_Buffer*min(BLER_Max);

Min_Num_Of_Errors_Sufficient_Statistics = 100;

file_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Files');
plot_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Plots');
code_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Code');
addpath(genpath(code_path));

is_visible = 'off';

if(m_SNR_EbN0==1)
    BLER = inf(9,length(SNR_Vec));
    BER = inf(9,length(SNR_Vec));
elseif(m_SNR_EbN0==2)    
    BLER = inf(8,length(SNR_Vec));
    BER = inf(8,length(SNR_Vec));
end

SNR_Type = 'Es/N0';

file_name = sprintf('BER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),strrep(sprintf('%s',datestr(now)),':','_'));
plot_name = sprintf('BER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),strrep(sprintf('%s',datestr(now)),':','_'));
save(fullfile(file_path,file_name),'N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
%Try


%%
%----------------------------gArikan BPSK----------------------------%

if(m_SNR_EbN0==1)

    Constellation_Type = 'BPSK';
    m_BPSK = 1;
    R_Code_Bits_BPSK = R_Channel_Bits/m_BPSK; % [Information bits / Total bits]

    [BLER(1,:),BER(1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,'Eb/N0','gArikan',1,1,0,Constellation_Type,m_BPSK,zeros(0),zeros(0),'BPSK',N,R_Code_Bits_BPSK,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
    save(fullfile(file_path,file_name),'N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

end

%%
Constellation_Type = 'QAM';
m_4QAM = 2;
R_Code_Bits_4QAM = R_Channel_Bits/m_4QAM; % [Information bits / Total bits]

[BLER((m_SNR_EbN0==1)+1,:),BER((m_SNR_EbN0==1)+1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,Constellation_Type,m_4QAM,zeros(0),zeros(0),'SP',N,R_Code_Bits_4QAM,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
Constellation_Type = 'QAM';
m_16QAM = 4;
R_Code_Bits_16QAM = R_Channel_Bits/m_16QAM; % [Information bits / Total bits]

[BLER((m_SNR_EbN0==1)+2,:),BER((m_SNR_EbN0==1)+2,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,Constellation_Type,m_16QAM,zeros(0),zeros(0),'SP',N,R_Code_Bits_16QAM,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
Constellation_Type = 'QAM';
m_256QAM = 8;
R_Code_Bits_256QAM = R_Channel_Bits/m_256QAM; % [Information bits / Total bits]

[BLER((m_SNR_EbN0==1)+3,:),BER((m_SNR_EbN0==1)+3,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,Constellation_Type,m_256QAM,zeros(0),zeros(0),'SP',N,R_Code_Bits_256QAM,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
Constellation_Type = 'PAM';
m_4PAM = 2;
R_Code_Bits_4PAM = R_Channel_Bits/m_4PAM; % [Information bits / Total bits]

[BLER((m_SNR_EbN0==1)+4,:),BER((m_SNR_EbN0==1)+4,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,Constellation_Type,m_4PAM,zeros(0),zeros(0),'SP',N,R_Code_Bits_4PAM,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
Constellation_Type = 'PAM';
m_16PAM = 4;
R_Code_Bits_16PAM = R_Channel_Bits/m_16PAM; % [Information bits / Total bits]

[BLER((m_SNR_EbN0==1)+5,:),BER((m_SNR_EbN0==1)+5,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,Constellation_Type,m_16PAM,zeros(0),zeros(0),'SP',N,R_Code_Bits_16PAM,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
Constellation_Type = 'PAM';
m_256PAM = 8;
R_Code_Bits_256PAM = R_Channel_Bits/m_256PAM; % [Information bits / Total bits]

[BLER((m_SNR_EbN0==1)+6,:),BER((m_SNR_EbN0==1)+6,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,Constellation_Type,m_256PAM,zeros(0),zeros(0),'SP',N,R_Code_Bits_256PAM,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
Constellation_Type = 'PSK';
m_16PSK = 4;
R_Code_Bits_16PSK = R_Channel_Bits/m_16PSK; % [Information bits / Total bits]

[BLER((m_SNR_EbN0==1)+7,:),BER((m_SNR_EbN0==1)+7,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,Constellation_Type,m_16PSK,zeros(0),zeros(0),'SP',N,R_Code_Bits_16PSK,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
Constellation_Type = 'PSK';
m_256PSK = 8;
R_Code_Bits_256PSK = R_Channel_Bits/m_256PSK; % [Information bits / Total bits]

[BLER((m_SNR_EbN0==1)+8,:),BER((m_SNR_EbN0==1)+8,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,'MLPC',1,1,0,Constellation_Type,m_256PSK,zeros(0),zeros(0),'SP',N,R_Code_Bits_256PSK,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
%----------------------------plot BLER----------------------------%

if(m_SNR_EbN0==1)
    legend_name = {sprintf('BPSK gArikan R=%s[information bits per channel use]',sprintf('%.3f',R_Code_Bits_BPSK)),sprintf('4QAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_4QAM)),sprintf('16QAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_16QAM)),sprintf('256QAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_256QAM)),sprintf('4PAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_4PAM)),sprintf('16PAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_16PAM)),sprintf('256PAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_256PAM)),sprintf('16PSK gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_16PSK)),sprintf('256PSK gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_256PSK))};
elseif(m_SNR_EbN0==2)
    legend_name = {sprintf('4QAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_4QAM)),sprintf('16QAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_16QAM)),sprintf('256QAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_256QAM)),sprintf('4PAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_4PAM)),sprintf('16PAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_16PAM)),sprintf('256PAM gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_256PAM)),sprintf('16PSK gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_16PSK)),sprintf('256PSK gArikan R=%s[bits per channel use]',sprintf('%.3f',R_Code_Bits_256PSK))};
end

title_name = sprintf('N=%d[bits] R=%s[information bits per channel use]',N,sprintf('%.3f',R_Channel_Bits));

save(fullfile(file_path,file_name),'plot_name','legend_name','title_name','N','R_Channel_Bits','m_SNR_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

Plot_BLER_or_BER_Comparison(false,SNR_Type,SNR_Vec_dB_EsN0,BER,length(BER(:,1)),Min_BLER,title_name,legend_name,plot_path,plot_name);

%%
sprintf('finish') 