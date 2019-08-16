close all;
%clear all;

%%
%----------------------------Parallel Computing----------------------------%

myCluster = parcluster('local');
NumOfWorkers = myCluster.NumWorkers;
try
   parpool local
end

%%
%--------------------------------N = 2^9--------------------------------%
%N = 2^9; % Number of bits in block

%R_Code = 7/8;
%SNR_Vec_dB = [?];

%R_Code = 3/4;
%SNR_Vec_dB = [?];

%R_Code = 1/2;
%SNR_Vec_dB = [?];

%R_Code = 1/4;
%SNR_Vec_dB = [?];

%R_Code = 1/8;
%SNR_Vec_dB = [?];

%%
%--------------------------------N = 2^10--------------------------------%
%N = 2^10; % Number of bits in block

%Compound Polar Codes 16QAM
%R_Code = 1/2;
%SNR_Vec_dB_EbN0 = [?];

%%
%--------------------------------N = 2^11--------------------------------%
% N = 2^11; % Number of bits in block
% 
% %List Decoding of Polar Codes BPSK
% R_Channel_Bits = 1/2; % [bits/channel use]
% SNR_Vec_dB_EbN0 = [1.75,2,2.25,2.5,2.75,3];

%%
%Try
N = 2^12; % Number of bits in block

SNR_Vec_dB_EbN0 = [0:0.5:3];

%%

Constellation_Type = 'BPSK';
Decode_Type = 'gArikan';
%Decode_Type = 'MLPC';
Lable_Type = 'BPSK';
%Lable_Type = 'SP';
m = 1;
R = 1/2;
R_Channel_Bits = m*R; % [bits/channel use]
if(m==1)
    m_SCL = 1;
else
    m_SCL = m./(2.^[log2(m):-1:0]);
end
R_Code_Bits = R_Channel_Bits/m; % [Information bits / Total bits]

%%

L_List = [1,2,4,8,16,32];

L_Optimize_List = ones(1,size(L_List,2));

GA_CRC_Length_List = zeros(1,size(L_List,2));
% GA_CRC_Length_List = 16*ones(1,size(L_List,2)); %crc16

SNR_Vec_dB_EsN0 = SNR_Vec_dB_EbN0 + 10*log10(R_Channel_Bits);
SNR_Vec = 10.^(SNR_Vec_dB_EsN0./10);
BLER_Max = 10^(-2);

BLER_Axis_Buffer = 10^-1;
Min_BLER = BLER_Axis_Buffer*min(BLER_Max);

Min_Num_Of_Errors_Sufficient_Statistics = 100;

file_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Files');
plot_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Plots');
code_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Code');
addpath(genpath(code_path));

is_visible = 'off';

BLER = inf(1+size(L_List,2)*size(m_SCL,2),size(SNR_Vec,2));
BER = inf(1+size(L_List,2)*size(m_SCL,2),size(SNR_Vec,2));

SNR_Type = 'Eb/N0';

file_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),strrep(sprintf('%s',datestr(now)),':','_'));
plot_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),strrep(sprintf('%s',datestr(now)),':','_'));

% if exist(fullfile(file_path,file_name))
%     save(fullfile(file_path,file_name),'N','R_Channel_Bits','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','L_List');
% else
%     fcreate
% end
%%
%Try

[BLER,BER] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,'Eb/N0',Decode_Type,L_List,L_Optimize_List,GA_CRC_Length_List,Constellation_Type,m,zeros(0),zeros(0),Lable_Type,N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);

legend_name{size(L_List,2)} = [];

for i=1:1:size(L_List,2)
    
    if(GA_CRC_Length_List(i)>=0)
        legend_name{i} = sprintf('BPSK L=%d GA CRC %d',L_List(i),GA_CRC_Length_List(i));
    else
        legend_name{i} = sprintf('BPSK L=%d',L_List(i));
    end
    
end

title_name = sprintf('N=%d[bits] R=%s[information bits per channel use]',N,sprintf('%.3f',R_Channel_Bits));

% save(fullfile(file_path,file_name),'plot_name','legend_name','title_name','N','m','m_SCL','R_Channel_Bits','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','L_List');

%Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,size(L_List,2)/2,Min_BLER,title_name,legend_name,plot_path,plot_name);
Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,size(L_List,2),Min_BLER,title_name,legend_name,plot_path,plot_name);
Plot_BLER_or_BER_Comparison(false,SNR_Type,SNR_Vec_dB_EbN0,BER,size(L_List,2),Min_BLER,title_name,legend_name,plot_path,plot_name);

%%

if(m==1)

    Constellation_Type = 'BPSK';
    m = 1;
    R_Code_Bits = R_Channel_Bits/m; % [Information bits / Total bits]

    [BLER(1:size(L_List,2)+1,:),BER(1:size(L_List,2)+1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,'Eb/N0',Decode_Type,[1,L_List],[1,L_Optimize_List],[1,GA_CRC_Length_List],Constellation_Type,m,zeros(0),zeros(0),Lable_Type,N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
    save(fullfile(file_path,file_name),'N','m','R_Channel_Bits','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','L_List','L_Optimize_List');

else
    
    for i=1:1:size(m_SCL,2)
    
        if(i==1)
            [BLER((i-1)*size(L_List,2)+1:i*size(L_List,2)+1,:),BER((i-1)*size(L_List,2)+1:i*size(L_List,2)+1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,'Eb/N0',Decode_Type,[1,L_List],[1,L_Optimize_List],[1,GA_CRC_Length_List],Constellation_Type,m,zeros(0),m_SCL(i),Lable_Type,N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
        else
            [BLER((i-1)*size(L_List,2)+2:i*size(L_List,2)+1,:),BER((i-1)*size(L_List,2)+2:i*size(L_List,2)+1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,'Eb/N0',Decode_Type,L_List,L_Optimize_List,GA_CRC_Length_List,Constellation_Type,m,zeros(0),m_SCL(i),Lable_Type,N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
        end
        save(fullfile(file_path,file_name),'N','m','m_SCL','R_Channel_Bits','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','L_List','L_Optimize_List');
    end
    
end
    
%%
%----------------------------plot BLER----------------------------%

legend_name{size(L_List,2)*size(m_SCL,2)+1} = [];

legend_name{1} = sprintf('%d%s',2^m,Constellation_Type);

for i=1:1:size(m_SCL,2)
    for j=1:1:size(L_List,2)
        if(m==1)
            if(L_List(j)==1)
                legend_name{(i-1)*size(L_List,2)+j+1} = sprintf('%d%s',2^m,Constellation_Type);
            else
                legend_name{(i-1)*size(L_List,2)+j+1} = sprintf('%d%s L=%d GA CRC %d',2^m,Constellation_Type,L_List(j),GA_CRC_Length_List(j));
            end
        else
            if(L_List(j)==1)
                legend_name{(i-1)*size(L_List,2)+j+1} = sprintf('%d%s %s %s',2^m,Constellation_Type,Decode_Type,Lable_Type);
            else
                legend_name{(i-1)*size(L_List,2)+j+1} = sprintf('%d%s %s %s L=%d(%d/%d) GA CRC %d',2^m,Constellation_Type,Decode_Type,Lable_Type,L_List(j),m_SCL(i),m,GA_CRC_Length_List(j));
            end
        end
    end
end

title_name = sprintf('N=%d[bits] R=%s[information bits per channel use]',N,sprintf('%.3f',R_Channel_Bits));

save(fullfile(file_path,file_name),'plot_name','legend_name','title_name','N','m','m_SCL','R_Channel_Bits','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','L_List');

Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,size(BLER,1),Min_BLER,title_name,legend_name,plot_path,plot_name);

%%
sprintf('finish') 