close all;
% clear all;

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
% N = 2*4^5; % Number of bits in block
% N = 2*4^4; % Number of bits in block

m = 4;
% m = 8;

% Constellation_Type = '16APSK-2-Rings';
% Constellation_Type = '16APSK-8-8';
% R_Channel_Bits = m*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [4,5,6];

% Constellation_Type = 'QAM';
% R_Channel_Bits = m*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [9.5,10,10.5,11,11.5,12,12.5];
% SNR_Vec_dB_EbN0 = [9.5,10,10.5];
% SNR_Vec_dB_EbN0 = [8];
% SNR_Vec_dB_EbN0 = [11];

% Constellation_Type = 'PAM';
% R_Channel_Bits = m*(3/4); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [13,13.5,14,14.5];
% SNR_Vec_dB_EbN0 = [14.5];

% Constellation_Type = 'QAM';
% R_Channel_Bits = m*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [6,6.5,7,7.5];

Constellation_Type = 'PAM';
R_Channel_Bits = m*(1/2); % [bits/channel use]
% SNR_Vec_dB_EbN0 = [9,10,11];
% SNR_Vec_dB_EbN0 = [9,9.5,10];
% SNR_Vec_dB_EbN0 = [8.5,9];
% SNR_Vec_dB_EbN0 = [10];
SNR_Vec_dB_EbN0 = [11];

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

% Decode_Type = 'MLPC';
% Lable_Type = 'SP';

% Decode_Type = 'BIMLPCM';
% Lable_Type = 'SP-Compound';

% Decode_Type = 'BIMLPCM_gRS4';
% Decode_Type = 'BIMLPCM_gArikan_GF4';
% Lable_Type = 'SP-Gray-GF4';

Decode_Type = 'BIPCM';
Lable_Type = 'Gray';

% Decode_Type = 'BIPCM_gRS4';
% Decode_Type = 'BIPCM_gArikan_GF4';
% Lable_Type = 'Gray-GF4';

m_BIPCM = m;
% m_BIPCM = m/2;
% m_BIPCM = 2;
% m_BIPCM = 1;
% m_BIPCM = 0;

is_sorted_by_N_m = true;
% is_sorted_by_N_m = false;

% is_sorted_by_Union_Bound = true;
is_sorted_by_Union_Bound = false;

if(and(is_sorted_by_N_m,is_sorted_by_Union_Bound))
    Num_Of_Labels_To_Compare = 7;
elseif(and(is_sorted_by_N_m,~is_sorted_by_Union_Bound))
    Num_Of_Labels_To_Compare = 5;
elseif(and(~is_sorted_by_N_m,is_sorted_by_Union_Bound))
    Num_Of_Labels_To_Compare = 3;
end

BLER = inf(Num_Of_Labels_To_Compare,length(SNR_Vec));
BER = inf(Num_Of_Labels_To_Compare,length(SNR_Vec));

SNR_Type = 'Eb/N0';

BLER_Constellation_Polar_Code = 10^(-2);

SNR_Constellation_Polar_Code_Vec = SNR_Vec(end);
SNR_Constellation_Union_Bound_Vec = SNR_Vec(end);

date_string = strrep(sprintf('%s',datestr(now)),':','_');
file_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%d%s_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),2^m,Constellation_Type,date_string);
plot_name = sprintf('BLER_N%d[bits]_R%s[information_bits_per_channel_use]_gArikan_%d%s_%s',N,strrep(sprintf('%.3f',R_Channel_Bits),'.','_'),2^m,Constellation_Type,date_string);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER');

%%
%Try

%load('D:\פאר\אוניברסיטה\תואר שני\תיזה - polar codes\MATLAB\Results\BIPCM Gray\16PAM R=0.75\Files\Same_Type_Labels_N=1024[bits]_R=0_750[bits_per_channel_use]_BIPCM_16PAM_All-Gray_gArikan.mat');
%load('D:\פאר\אוניברסיטה\תואר שני\תיזה - polar codes\MATLAB\Results\BIPCM Gray\16QAM R=0.75\Files\Same_Type_Labels_N=1024[bits]_R=0_750[bits_per_channel_use]_BIPCM_16QAM_All-Gray_gArikan.mat');

%%

[Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound,Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound,Sorted_By_Union_Bound] = Sort_Labelings(NumOfWorkers,Decode_Type,Constellation_Type,m,m_BIPCM,sprintf('All-%s',Lable_Type),N,R_Code_Bits,[SNR_Constellation_Polar_Code_Vec,SNR_Constellation_Union_Bound_Vec],BLER_Constellation_Polar_Code,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible,is_sorted_by_N_m,is_sorted_by_Union_Bound);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound','Sorted_By_Union_Bound');

% [BLER(1,:),BER(1,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,Decode_Type,1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),Lable_Type,N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
% save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound','Sorted_By_Union_Bound');

[BLER(2,:),BER(2,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,Decode_Type,1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound(1,:),N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound','Sorted_By_Union_Bound');

[BLER(3,:),BER(3,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,Decode_Type,1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound(end,:),N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound','Sorted_By_Union_Bound');

[BLER(4,:),BER(4,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,Decode_Type,1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound(1,:),N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound','Sorted_By_Union_Bound');

[BLER(5,:),BER(5,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,Decode_Type,1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound(end,:),N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound','Sorted_By_Union_Bound');

[BLER(6,:),BER(6,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,Decode_Type,1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),Sorted_By_Union_Bound(1,:),N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound','Sorted_By_Union_Bound');

[BLER(7,:),BER(7,:)] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,Decode_Type,1,1,0,Constellation_Type,m,m_BIPCM,zeros(0),Sorted_By_Union_Bound(end,:),N,R_Code_Bits,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible);
save(fullfile(file_path,file_name),'N','R_Code_Bits','Constellation_Type','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound','Sorted_By_Union_Bound');

%%
%----------------------------plot BLER----------------------------%

legend_name{Num_Of_Labels_To_Compare} = [];

legend_name{1} = sprintf('%s %d%s %s regular',Decode_Type,2^m,Constellation_Type,Lable_Type);

if(and(is_sorted_by_N_m,is_sorted_by_Union_Bound))

    legend_name{2} = sprintf('%s %d%s %s best N=m PDF lower bound',Decode_Type,2^m,Constellation_Type,Lable_Type);

    legend_name{3} = sprintf('%s %d%s %s worst N=m PDF lower bound',Decode_Type,2^m,Constellation_Type,Lable_Type);

    legend_name{4} = sprintf('%s %d%s %s best N=m PDF upper bound',Decode_Type,2^m,Constellation_Type,Lable_Type);

    legend_name{5} = sprintf('%s %d%s %s worst N=m PDF upper bound',Decode_Type,2^m,Constellation_Type,Lable_Type);
    
    legend_name{6} = sprintf('%s %d%s %s best union bound',Decode_Type,2^m,Constellation_Type,Lable_Type);

    legend_name{7} = sprintf('%s %d%s %s worst union bound',Decode_Type,2^m,Constellation_Type,Lable_Type);
    
elseif(and(is_sorted_by_N_m,~is_sorted_by_Union_Bound))
    
    legend_name{2} = sprintf('%s %d%s %s best N=m PDF lower bound',Decode_Type,2^m,Constellation_Type,Lable_Type);

    legend_name{3} = sprintf('%s %d%s %s worst N=m PDF lower bound',Decode_Type,2^m,Constellation_Type,Lable_Type);

    legend_name{4} = sprintf('%s %d%s %s best N=m PDF upper bound',Decode_Type,2^m,Constellation_Type,Lable_Type);

    legend_name{5} = sprintf('%s %d%s %s worst N=m PDF upper bound',Decode_Type,2^m,Constellation_Type,Lable_Type);
    
elseif(and(~is_sorted_by_N_m,is_sorted_by_Union_Bound))
    
    legend_name{2} = sprintf('%s %d%s %s best union bound',Decode_Type,2^m,Constellation_Type,Lable_Type);

    legend_name{3} = sprintf('%s %d%s %s worst union bound',Decode_Type,2^m,Constellation_Type,Lable_Type);
    
end

    g_Arikan = [1,0;1,1];
   
    zero = 0;
    one = 1;
    alpha = 2;
    alpha2 = 3;
    
    alpha_x_alpha = alpha*alpha;
    alpha_x_alpha2 = alpha*alpha2;
    alpha2_x_alpha2 = alpha2*alpha2;
    
    g_RS4 = [one zero zero zero; one one zero zero; alpha2 alpha one zero; one one one alpha];
    g_Arikan_GF4 = [one zero; one alpha];

if(strcmp(Decode_Type,'BIPCM'))
    
    g0 = g_Arikan;
    for i=1:1:log2(m)-1
        g0 = kron(g0,g_Arikan);
    end
    g0 = bitrevorder(g0);
   
elseif(strcmp(Decode_Type,'BIPCM_gRS4'))
        
    g0 = g_RS4;
    for i=1:1:log2(m/2)/log2(size(g_RS4,1))-1
        g0 = kron(g0,g_RS4);
        g0(g0==alpha_x_alpha) = alpha2;
        g0(g0==alpha_x_alpha2) = one;
        g0(g0==alpha2_x_alpha2) = alpha;
    end
    g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,4),'RS4');

elseif(strcmp(Decode_Type,'BIPCM_gArikan_GF4'))
            
    g0 = g_Arikan_GF4;
    for i=1:1:log2(m/2)/log2(size(g_Arikan_GF4,1))-1
        g0 = kron(g0,g_Arikan_GF4);
        g0(g0==alpha_x_alpha) = alpha2;
        g0(g0==alpha_x_alpha2) = one;
        g0(g0==alpha2_x_alpha2) = alpha;
    end
    g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,2),'Arikan_GF4');
    
elseif(strcmp(Decode_Type,'MLPC'))
    
    g0 = 1;

elseif(strcmp(Decode_Type,'MLPC_gRS4'))

    g0 = 1;
    
elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
   
    g0 = 1;
    
elseif(strcmp(Decode_Type,'BIMLPCM'))

    g0 = g_Arikan;
    for i=1:1:log2(m)-1
        g0 = kron(g0,g_Arikan);
    end
    g0 = bitrevorder(g0);
    
elseif(strcmp(Decode_Type,'BIMLPCM_gRS4'))
       
    g0 = g_RS4;
    for i=1:1:log2(m)/log2(2*size(g_RS4,1))-1
        g0 = kron(g0,g_RS4);
        g0(g0==alpha_x_alpha) = alpha2;
        g0(g0==alpha_x_alpha2) = one;
        g0(g0==alpha2_x_alpha2) = alpha;
    end
    g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,4),'RS4');
                
elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4'))
            
    g0 = g_Arikan_GF4;
    for i=1:1:log2(m)/log2(2*size(g_Arikan_GF4,1))-1
        g0 = kron(g0,g_Arikan_GF4);
        g0(g0==alpha_x_alpha) = alpha2;
        g0(g0==alpha_x_alpha2) = one;
        g0(g0==alpha2_x_alpha2) = alpha;
    end
    g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,2),'Arikan_GF4');
            
end

title_name = sprintf('N=%d[bits] R=%s[information bits per channel use]',N,sprintf('%.3f',R_Channel_Bits));

save(fullfile(file_path,file_name),'plot_name','legend_name','title_name','N','R_Code_Bits','Constellation_Type','m','m_BIPCM','SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','BLER','BER','Sorted_By_Union_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound'); 

Plot_BLER_or_BER_Comparison(true,SNR_Type,SNR_Vec_dB_EbN0,BLER,length(BLER(:,1)),Min_BLER,title_name,legend_name,plot_path,plot_name);

Regular_Labeling = Make_Constellation(Constellation_Type,m,Lable_Type,g0);

if(and(is_sorted_by_N_m,is_sorted_by_Union_Bound))

    N_m_PDF(NumOfWorkers,Decode_Type,Constellation_Type,m,m_BIPCM,[Regular_Labeling;Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound(1,:);Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound(end,:);Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound(1,:);Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound(end,:);Sorted_By_Union_Bound(1,:);Sorted_By_Union_Bound(end,:)],Lable_Type,R_Code_Bits,SNR_Constellation_Polar_Code_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,legend_name,'on');
    
elseif(and(is_sorted_by_N_m,~is_sorted_by_Union_Bound))
    
    N_m_PDF(NumOfWorkers,Decode_Type,Constellation_Type,m,m_BIPCM,[Regular_Labeling;Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound(1,:);Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound(end,:);Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound(1,:);Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound(end,:)],Lable_Type,R_Code_Bits,SNR_Constellation_Polar_Code_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,legend_name,'on');

elseif(and(~is_sorted_by_N_m,is_sorted_by_Union_Bound))
    
    N_m_PDF(NumOfWorkers,Decode_Type,Constellation_Type,m,m_BIPCM,[Regular_Labeling;Sorted_By_Union_Bound(1,:);Sorted_By_Union_Bound(end,:)],Lable_Type,R_Code_Bits,SNR_Constellation_Polar_Code_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,legend_name,'on');
    
end

%%
sprintf('finish')
    