% close all;
clear all;

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


%% list decoder related parameters
L_List = 1;% [1,2,4,8,16,32];
CRC_length = 0;
L_Optimize_List = ones(1,size(L_List,2));
GA_CRC_Length_List = CRC_length*ones(1,size(L_List,2)); %crc16

SNR_Vec_dB_EbN0 = 1.5:0.5:3.5;

%% code related parameters
N = 2^10; % Number of bits in block
Constellation_Type = 'BPSK';
Decode_Type = 'gArikan';
Lable_Type = 'BPSK';
m = 1;
R = 0.5; % desired code rate
Rtot = (R*N - CRC_length)/N; % total code rate after CRC addition
R_Channel_Bits = m*R; % [bits/channel use]
m_SCL = 1;
R_Code_Bits = R_Channel_Bits/m; % [Information bits / Total bits]


%% Polar Code related parameters
U_N = N;
Is_GF4 = false;
g_Arikan = [1,0;1,1];
M = 2^m;

% Generate polar encoding matrix
g = g_Arikan;
for i=1:1:log2(N)-1
    g = kron(g,g_Arikan);
end
g = bitrevorder(g);

g0 = 1;

Constellation_Mapping_Array = [-1,1];

%% Arrays initialization
SNR_Vec_dB_EsN0 = SNR_Vec_dB_EbN0 + 10*log10(Rtot);
SNR_Vec = 10.^(SNR_Vec_dB_EsN0./10);
BLER_Max = 10^(-2);

BLER_Axis_Buffer = 10^-1;
Min_BLER = BLER_Axis_Buffer*min(BLER_Max);

Min_Num_Of_Errors_Sufficient_Statistics = 100;

file_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Files');
plot_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'Plots');
code_path = fullfile(fileparts(mfilename('fullpath')),'PeerCode');
addpath(genpath(code_path));

is_visible = 'off';

BLER = inf(1+size(L_List,2)*size(m_SCL,2),size(SNR_Vec,2));
BER = inf(1+size(L_List,2)*size(m_SCL,2),size(SNR_Vec,2));

SNR_Type = 'Eb/N0';

file_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),strrep(sprintf('%s',datestr(now)),':','_'));
plot_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),strrep(sprintf('%s',datestr(now)),':','_'));
title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));

%%
% actual simulation starts here
SNR_Vec_dB_EsN0 = 10.*log10(SNR_Vec);
SNR_Vec_dB_EbN0 = SNR_Vec_dB_EsN0 - 10*log10(m) - 10*log10(R);

SER_Vec = zeros(size(SNR_Vec,2),U_N);
Is_Frozen_Bit_Index_Vec = zeros(size(SNR_Vec,2),N);
Is_Frozen_Bit_Index_Vec = Is_Frozen_Bit_Index_Vec == ones(size(SNR_Vec,2),N);
Is_Frozen_Bit_Index_List_Vec = zeros(size(L_List,2),size(SNR_Vec,2),N);
Is_Frozen_Bit_Index_List_Diff_Vec = zeros(size(L_List,2),size(SNR_Vec,2));
Is_Frozen_Bit_Index_temp_Vec = zeros(1,N);
Genie_Aided_LLR_Avg_Vec = zeros(size(SNR_Vec,2),U_N);
Genie_Aided_LLR_Var_Vec = zeros(size(SNR_Vec,2),U_N);
Estimated_SER_Vec = zeros(size(SNR_Vec,2),U_N);
BLER = zeros(size(L_List,2),size(SNR_Vec,2));
BER = zeros(size(L_List,2),size(SNR_Vec,2));
First_Symbol = zeros(size(L_List,2),size(SNR_Vec,2),U_N);
temp_First_Symbol = zeros(1,U_N);
temp_LLR_Avg = zeros(1,U_N);
temp_LLR_Var = zeros(1,U_N);
LLR_Avg = zeros(size(L_List,2),size(SNR_Vec,2),U_N);
LLR_Var = zeros(size(L_List,2),size(SNR_Vec,2),U_N);
Is_Change_Frozen_Bits_Index_Vec = zeros(size(SNR_Vec,2)-1,N);
Estimated_Bhattacharyya = zeros(size(SNR_Vec,2),U_N);
is_high_SNR = zeros(1,size(L_List,2));

for i=1:1:size(SNR_Vec,2)
    
    if(sum(is_high_SNR)==size(is_high_SNR,2))
        break
    end
    
    %%  Estimate the best set of frozen bits
    [Estimated_Bhattacharyya(i,:),SER_Vec(i,:),Is_Frozen_Bit_Index_Vec(i,:),~,Genie_Aided_LLR_Avg_Vec(i,:),Genie_Aided_LLR_Var_Vec(i,:),...
        Estimated_SER_Vec(i,:),~,~,~,~,~] = listDecoderSim(NumOfWorkers,true,g,g0,Decode_Type,1,0,Constellation_Type,m,[],m_SCL,Constellation_Mapping_Array,N,R,BLER_Max,...
        SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_Vec(i,:));
    
    %     Plot_Index_SER(m,Estimated_Bhattacharyya(i,:),SER_Vec(i,:),Estimated_SER_Vec(i,:),Genie_Aided_LLR_Avg_Vec(i,:),sqrt(Genie_Aided_LLR_Var_Vec(i,:)),...
    %         Is_Frozen_Bit_Index_Vec(i,:),title_name,'BPSK gArikan',plot_path,sprintf('%s SNR=%s',plot_name,strrep(sprintf('%.3f',SNR_Vec_dB_EbN0(1,i)),'.','_')),is_visible);
    %
    
    Is_Frozen_Bit_Index_Vec(i,:) = Is_Frozen_Bit_Index_Vec(i,:) == ones(1,N);
    
    for j=1:1:size(L_List,2)
        
        if(is_high_SNR(1,j))
            continue
        end
        
        if(L_List(j)==1)
            
            Is_Frozen_Bit_Index_List_Vec(j,i,:) = Is_Frozen_Bit_Index_Vec(i,:);
            
            Is_Frozen_Bit_Index_temp_Vec(1:N) = Is_Frozen_Bit_Index_List_Vec(j,i,:);
            %%  BER Simulation with the calculated frozen bits set
            
            [~,~,~,~,~,~,~,BLER(1,i),BER(1,i),First_Symbol(1,i,:),LLR_Avg(1,i,:),LLR_Var(1,i,:)] = listDecoderSim(NumOfWorkers,false,g,g0,Decode_Type,L_List(j),GA_CRC_Length_List(j),Constellation_Type,m,[],...
                m_SCL,Constellation_Mapping_Array,N,R,BLER_Max,SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_temp_Vec);
            
            if(i>1)
                Is_Change_Frozen_Bits_Index_Vec(i-1,:) = Is_Frozen_Bit_Index_Vec(i-1,:)~=Is_Frozen_Bit_Index_Vec(i,:);
                sprintf('%d frozen bits has changed',sum(Is_Change_Frozen_Bits_Index_Vec(i-1,:))/2)
            end
            
        else
            
            if(L_Optimize_List(j)==1)
                
                Is_Frozen_Bit_Index_List_Vec(j,i,:) = Is_Frozen_Bit_Index_Vec(i,:);
                
            else
                
                [Sorted_SER_Vec,Sorted_SER_Index_Vec] = sort(wrev(SER_Vec(i,:)));
                Sorted_SER_Index_Vec = size(SER_Vec,2) - Sorted_SER_Index_Vec + 1;
                
                [Is_Frozen_Bit_Index_List_Vec(j,i,:),Is_Frozen_Bit_Index_List_Diff_Vec(j,i)] = Find_Frozen_Bits_List(Sorted_SER_Vec,Sorted_SER_Index_Vec,R,L_Optimize_List(j));
            end
            
            Is_Frozen_Bit_Index_temp_Vec(1:N) = Is_Frozen_Bit_Index_List_Vec(j,i,:);
            
            [~,~,~,~,~,~,~,BLER(j,i),BER(j,i),First_Symbol(j,i,:),LLR_Avg(j,i,:),LLR_Var(j,i,:)] = listDecoderSim(NumOfWorkers,false,g,g0,Decode_Type,L_List(j),GA_CRC_Length_List(j),Constellation_Type,m,[],m_SCL,Constellation_Mapping_Array,N,R,BLER_Max,SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_temp_Vec);
            
        end
        
        temp_First_Symbol(1,:) = First_Symbol(j,i,:);
        temp_LLR_Avg(1,:) = LLR_Avg(j,i,:);
        temp_LLR_Var(1,:) = LLR_Var(j,i,:);
        
        if(L_List(j)~=1)
            if(GA_CRC_Length_List(j)>=0)
                if(or(m_SCL~=m,m==1))
                    plot_name_temp = sprintf('L=%d_GA_CRC_%d_%s',L_List(j),GA_CRC_Length_List(j),plot_name);
                    legend_name_temp = sprintf('L=%d(%d/%d) with GA CRC %d %s',L_List(j),m_SCL,m,GA_CRC_Length_List(j),'gArikan');
                else
                    plot_name_temp = sprintf('L=%d(%dof%d)_GA_CRC_%d_%s',L_List(j),m_SCL,m,GA_CRC_Length_List(j),plot_name);
                    legend_name_temp = sprintf('L=%d with GA CRC %d %s',L_List(j),GA_CRC_Length_List(j),'gArikan');
                end
            else
                if(or(m_SCL~=m,m==1))
                    plot_name_temp = sprintf('L=%d_%s',L_List(j),plot_name);
                    legend_name_temp = sprintf('L=%d(%d/%d) %s',L_List(j),m_SCL,m,'gArikan');
                else
                    plot_name_temp = sprintf('L=%d(%dof%d)_%s',L_List(j),m_SCL,m,plot_name);
                    legend_name_temp = sprintf('L=%d %s',L_List(j),'gArikan');
                end
            end
        else
            plot_name_temp = plot_name;
            legend_name_temp = strcat(strcat('gArikan L = ',num2str(L_List(j))));
        end
        
        %         Plot_Simulation_Results(m,SNR_Type,SNR_Vec_dB_EbN0(1,1:i),BLER(j,1:i),Min_BLER,BER(j,1:i),temp_First_Symbol,temp_LLR_Avg,sqrt(temp_LLR_Var),Is_Frozen_Bit_Index_Vec(i,:),title_name,legend_name_temp,plot_path,sprintf('%s SNR=%s',plot_name_temp,strrep(sprintf('%.3f',SNR_Vec_dB_EbN0(1,i)),'.','_')),is_visible);
        
        try
            save(fullfile(file_path,file_name),'SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','SER_Vec','Genie_Aided_LLR_Avg_Vec','Genie_Aided_LLR_Var_Vec','Estimated_SER_Vec','BLER','BER','First_Symbol','LLR_Avg','LLR_Var','Is_Frozen_Bit_Index_Vec','Is_Change_Frozen_Bits_Index_Vec','Is_Frozen_Bit_Index_List_Vec','Is_Frozen_Bit_Index_List_Diff_Vec');
        end
        
        if(and(BLER(j,i)<BLER_Max,j<size(L_List,2)))
            is_high_SNR(1,j:end) = 1;
            break;
        end
        
        if(and(BLER(j,i)<BLER_Max*10,i<size(SNR_Vec,2)))
            is_high_SNR(1,j) = 1;
        end
        
    end
    
    try
        save(fullfile(file_path,file_name),'SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','SER_Vec','Genie_Aided_LLR_Avg_Vec','Genie_Aided_LLR_Var_Vec','Estimated_SER_Vec','BLER','BER','First_Symbol','LLR_Avg','LLR_Var','Is_Frozen_Bit_Index_Vec','Is_Change_Frozen_Bits_Index_Vec','Is_Frozen_Bit_Index_List_Vec','Is_Frozen_Bit_Index_List_Diff_Vec');
    end
    
end

sprintf('%d frozen bits has changed\n',sum(Is_Change_Frozen_Bits_Index_Vec,2)./2)

title_name = sprintf('N=%d[bits] R=%s[information bits per channel use]',N,sprintf('%.3f',R_Channel_Bits));

% save(fullfile(file_path,file_name),'plot_name','legend_name','title_name','N','m','m_SCL','R_Channel_Bits','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','L_List');

%Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,size(L_List,2)/2,Min_BLER,title_name,legend_name,plot_path,plot_name);
% Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,size(L_List,2),Min_BLER,title_name,'gArikan',plot_path,plot_name);
% Plot_BLER_or_BER_Comparison(false,SNR_Type,SNR_Vec_dB_EbN0,BER,size(L_List,2),Min_BLER,title_name,'gArikan',plot_path,plot_name);


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

% save(fullfile(file_path,file_name),'plot_name','legend_name','title_name','N','m','m_SCL','R_Channel_Bits','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','L_List');

Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,size(BLER,1),Min_BLER,title_name,legend_name,plot_path,plot_name);

%%
sprintf('finish')