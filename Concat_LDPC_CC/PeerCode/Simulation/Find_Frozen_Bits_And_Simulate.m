function [BLER,BER] = Find_Frozen_Bits_And_Simulate(NumOfWorkers,SNR_Type,Decode_Type,L_List,L_Optimize_List,GA_CRC_Length_List,Constellation_Type,m,m_BIPCM,m_SCL,Lable_Type,N,R,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible)

    g_Arikan = [1,0;1,1];
    
    zero = 0;
    one = 1;
    alpha = 2;
    alpha2 = 3;
    
    alpha_x_alpha = alpha*alpha;
    alpha_x_alpha2 = alpha*alpha2;
    alpha2_x_alpha2 = alpha2*alpha2;
    
%     g2_Arikan = gf([one,zero;one,alpha],2);
    
%     g_RS4 = gf([one zero zero zero; one one zero zero; alpha2 alpha one zero; one one one alpha],2);
    g_RS4 = [one zero zero zero; one one zero zero; alpha2 alpha one zero; one one one alpha];

    g_Arikan_GF4 = [one zero; one alpha];
    
    M = 2^m;
    
    if(strcmp(Decode_Type,'gRS4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'MLPC_gRS4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIPCM_gRS4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIMLPCM_gRS4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'gArikan_GF4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIPCM_gArikan_GF4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4'))
        U_N = N/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'MLPC_gRS4_gUVV4'))
        U_N = N/2;
        Is_GF4 = true;
    else %gArikan
        U_N = N;
        Is_GF4 = false;
    end
    
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

    if(R>1)
        BLER = ones(size(L_List,2),size(SNR_Vec,2));
        BER = ones(size(L_List,2),size(SNR_Vec,2));
        return;
    end

    if(M==2)
        
        if(strcmp(Decode_Type,'gArikan'))
        
            g = g_Arikan;
            for i=1:1:log2(N)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);

            g0 = 1;

            Constellation_Mapping_Array = [-1,1];

            plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BPSK_gArikan',N,strrep(sprintf('%.3f',R),'.','_'));
            file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BPSK_gArikan',N,strrep(sprintf('%.3f',R),'.','_'));
            legend_name = sprintf('BPSK gArikan');
            title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
        elseif(strcmp(Decode_Type,'gRS4'))
            
            g = g_RS4;
            for i=1:1:log2(U_N)/log2(size(g_RS4,1))-1
                g = kron(g,g_RS4);
                g(g==alpha_x_alpha) = alpha2;
                g(g==alpha_x_alpha2) = one;
                g(g==alpha2_x_alpha2) = alpha;
            end
            g = Convert_Matrix_GF4_to_Matrix(digitrevorder(g,4),'RS4');
            
            g0 = 1;
            
            Constellation_Mapping_Array = [-1,1];

            plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BPSK_gRS4',N,strrep(sprintf('%.3f',R),'.','_'));
            file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BPSK_gRS4',N,strrep(sprintf('%.3f',R),'.','_'));
            legend_name = sprintf('BPSK gRS4');
            title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
        elseif(strcmp(Decode_Type,'gArikan_GF4'))
            
            g = g_Arikan_GF4;
            for i=1:1:log2(U_N)/log2(size(g_Arikan_GF4,1))-1
                g = kron(g,g_Arikan_GF4);
                g(g==alpha_x_alpha) = alpha2;
                g(g==alpha_x_alpha2) = one;
                g(g==alpha2_x_alpha2) = alpha;
            end
            g = Convert_Matrix_GF4_to_Matrix(digitrevorder(g,2),'Arikan_GF4');
            
            g0 = 1;
            
            Constellation_Mapping_Array = [-1,1];
            
            plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BPSK_gArikan_GF4',N,strrep(sprintf('%.3f',R),'.','_'));
            file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BPSK_gArikan_GF4',N,strrep(sprintf('%.3f',R),'.','_'));
            legend_name = sprintf('BPSK gArikan GF4');
            title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
            
        end
        
    elseif(mod(log(M)/log(4),1)==0)
        
        if(strcmp(Decode_Type,'Multi_Dimensional_MLPC'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = 1;
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),Constellation_Type);
                legend_name = sprintf('MLPC %s gArikan',Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),Constellation_Type,Lable_Type);
                legend_name = sprintf('MLPC %s %s gArikan',Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            end
            
        elseif(strcmp(Decode_Type,'Multi_Dimensional_BIPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
                        
            g0 = g_Arikan;
            for i=1:1:log2(m)-1
                g0 = kron(g0,g_Arikan);
            end
            g0 = bitrevorder(g0);
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),Constellation_Type);
                legend_name = sprintf('BIPCM %s gArikan',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),Constellation_Type,Lable_Type);
                legend_name = sprintf('BIPCM %s %s gArikan',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
            end
            
        elseif(strcmp(Decode_Type,'Multi_Dimensional_BIMLPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = g_Arikan;
            for i=1:1:log2(m)-1
                g0 = kron(g0,g_Arikan);
            end
            g0 = bitrevorder(g0);
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('BIMLPCM %d%s gArikan',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('BIMLPCM %d%s %s gArikan',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
            end
            
        elseif(strcmp(Decode_Type,'MLPC'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = 1;
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('MLPC %d%s gArikan',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('MLPC %d%s %s gArikan',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            end
        
        elseif(strcmp(Decode_Type,'MLPC_gRS4'))
            
            g = g_RS4;
            for i=1:1:log2(N/m)/log2(size(g_RS4,1))-1
                g = kron(g,g_RS4);
                g(g==alpha_x_alpha) = alpha2;
                g(g==alpha_x_alpha2) = one;
                g(g==alpha2_x_alpha2) = alpha;
            end
            g = Convert_Matrix_GF4_to_Matrix(digitrevorder(g,4),'RS4');
            
            g0 = 1;
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_gRS4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_gRS4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('MLPC %d%s gRS4',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gRS4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gRS4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('MLPC %d%s %s gRS4',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            end
           
        elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
            
            g = g_Arikan_GF4;
            for i=1:1:log2(N/m)/log2(size(g_Arikan_GF4,1))-1
                g = kron(g,g_Arikan_GF4);
                g(g==alpha_x_alpha) = alpha2;
                g(g==alpha_x_alpha2) = one;
                g(g==alpha2_x_alpha2) = alpha;
            end
            g = Convert_Matrix_GF4_to_Matrix(digitrevorder(g,2),'Arikan_GF4');
            
            g0 = 1;
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_gArikan_GF4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_gArikan_GF4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('MLPC %d%s gArikan GF4',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gArikan_GF4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gArikan_GF4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('MLPC %d%s %s gArikan GF4',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            end
            
        elseif(strcmp(Decode_Type,'MLPC_gRS4_gUVV4'))
            
            g1 = g_RS4;
            for i=1:1:log2(N/m)/log2(size(g_RS4,1))-1
                g1 = kron(g1,g_RS4);
                g1(g1==alpha_x_alpha) = alpha2;
                g1(g1==alpha_x_alpha2) = one;
                g1(g1==alpha2_x_alpha2) = alpha;
            end
            
            g2 = g_Arikan_GF4;
            for i=1:1:log2(N/m)/log2(size(g_Arikan_GF4,1))-1
                g2 = kron(g2,g_Arikan_GF4);
                g2(g2==alpha_x_alpha) = alpha2;
                g2(g2==alpha_x_alpha2) = one;
                g2(g2==alpha2_x_alpha2) = alpha;
            end
            
            g = [Convert_Matrix_GF4_to_Matrix(digitrevorder(g1,2),'RS4');Convert_Matrix_GF4_to_Matrix(digitrevorder(g2,2),'Arikan_GF4')];
            
            g0 = 1;
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_gRS4_gUVV4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_gRS4_gUVV4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('MLPC %d%s gRS4-gUVV4',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gRS4_gUVV4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gRS4_gUVV4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('MLPC %d%s %s gRS4-gUVV4',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            end
            
        elseif(strcmp(Decode_Type,'Compound_MLPC'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = g_Arikan;
            for i=1:1:log2(m)-1
                g0 = kron(g0,g_Arikan);
            end
            g0 = bitrevorder(g0);
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_Compound_MLPC_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_Compound_MLPC_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('Compound-MLPC %d%s gArikan',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_Compound_MLPC_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_Compound_MLPC_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('Compound-MLPC %d%s %s gArikan',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            end
            
        elseif(strcmp(Decode_Type,'Separated_BIPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = 1;
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_Separated_BIPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_Separated_BIPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('Separated-BIPCM %d%s gArikan',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_Separated_BIPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_Separated_BIPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('Separated-BIPCM %d%s %s gArikan',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            end
            
        elseif(strcmp(Decode_Type,'BIPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
                        
            g0 = g_Arikan;
            for i=1:1:log2(m)-1
                g0 = kron(g0,g_Arikan);
            end
            g0 = bitrevorder(g0);
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('BIPCM %d%s gArikan',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('BIPCM %d%s %s gArikan',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
            end
        
        elseif(strcmp(Decode_Type,'BIPCM_gRS4'))
            
            g = g_RS4;
            for i=1:1:log2(U_N)/log2(size(g_RS4,1))-1
                g = kron(g,g_RS4);
                g(g==alpha_x_alpha) = alpha2;
                g(g==alpha_x_alpha2) = one;
                g(g==alpha2_x_alpha2) = alpha;
            end
            g = Convert_Matrix_GF4_to_Matrix(digitrevorder(g,4),'RS4');
                        
            g0 = g_RS4;
            for i=1:1:log2(m/2)/log2(size(g_RS4,1))-1
                g0 = kron(g0,g_RS4);
                g0(g0==alpha_x_alpha) = alpha2;
                g0(g0==alpha_x_alpha2) = one;
                g0(g0==alpha2_x_alpha2) = alpha;
            end
            g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,4),'RS4');
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_gRS4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_gRS4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('BIPCM %d%s gRS4',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_%s_gRS4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_%s_gRS4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('BIPCM %d%s %s gRS4',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
            end
           
        elseif(strcmp(Decode_Type,'BIPCM_gArikan_GF4'))
            
            g = g_Arikan_GF4;
            for i=1:1:log2(U_N)/log2(size(g_Arikan_GF4,1))-1
                g = kron(g,g_Arikan_GF4);
                g(g==alpha_x_alpha) = alpha2;
                g(g==alpha_x_alpha2) = one;
                g(g==alpha2_x_alpha2) = alpha;
            end
            g = Convert_Matrix_GF4_to_Matrix(digitrevorder(g,2),'Arikan_GF4');
                        
            g0 = g_Arikan_GF4;
            for i=1:1:log2(m/2)/log2(size(g_Arikan_GF4,1))-1
                g0 = kron(g0,g_Arikan_GF4);
                g0(g0==alpha_x_alpha) = alpha2;
                g0(g0==alpha_x_alpha2) = one;
                g0(g0==alpha2_x_alpha2) = alpha;
            end
            g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,2),'Arikan_GF4');
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_gArikan_GF4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_gArikan_GF4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('BIPCM %d%s gArikan GF4',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_%s_gArikan_GF4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_%s_gArikan_GF4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('BIPCM %d%s %s gArikan GF4',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
            end
            
        elseif(strcmp(Decode_Type,'Compound_BIPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = g_Arikan;
            g0 = bitrevorder(g0);
            g0 = kron(eye(m/2),g0);
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_Compound_BIPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_Compound_BIPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('Compound-BIPCM %d%s gArikan',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_Compound_BIPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_Compound_BIPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('Compound-BIPCM %d%s %s gArikan',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));        
          
            end
            
         elseif(strcmp(Decode_Type,'BIMLPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = g_Arikan;
            for i=1:1:log2(m)-1
                g0 = kron(g0,g_Arikan);
            end
            g0 = bitrevorder(g0);
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('BIMLPCM %d%s gArikan',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('BIMLPCM %d%s %s gArikan',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
            end
        
        elseif(strcmp(Decode_Type,'BIMLPCM_BIPCM_MLPC'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = g_Arikan;
            for i=1:1:log2(m)-1
                g0 = kron(g0,g_Arikan);
            end
            g0 = bitrevorder(g0);
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_gArikan',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('BIMLPCM %d%s gArikan',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gArikan',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('BIMLPCM %d%s %s gArikan',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
            end
            
        elseif(strcmp(Decode_Type,'BIMLPCM_gRS4'))
            
            g = g_RS4;
            for i=1:1:log2(N/m)/log2(size(g_RS4,1))-1
                g = kron(g,g_RS4);
                g(g==alpha_x_alpha) = alpha2;
                g(g==alpha_x_alpha2) = one;
                g(g==alpha2_x_alpha2) = alpha;
            end
            g = Convert_Matrix_GF4_to_Matrix(digitrevorder(g,4),'RS4');
            
            g0 = g_RS4;
            for i=1:1:log2(m/2)/log2(size(g_RS4,1))-1
                g0 = kron(g0,g_RS4);
                g0(g0==alpha_x_alpha) = alpha2;
                g0(g0==alpha_x_alpha2) = one;
                g0(g0==alpha2_x_alpha2) = alpha;
            end
            g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,4),'RS4');
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_gRS4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_gRS4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('BIMLPCM %d%s gRS4',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gRS4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gRS4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('BIMLPCM %d%s %s gRS4',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
            end
           
        elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4'))
            
            g = g_Arikan_GF4;
            for i=1:1:log2(N/m)/log2(size(g_Arikan_GF4,1))-1
                g = kron(g,g_Arikan_GF4);
                g(g==alpha_x_alpha) = alpha2;
                g(g==alpha_x_alpha2) = one;
                g(g==alpha2_x_alpha2) = alpha;
            end
            g = Convert_Matrix_GF4_to_Matrix(digitrevorder(g,2),'Arikan_GF4');
            
            g0 = g_Arikan_GF4;
            for i=1:1:log2(m/2)/log2(size(g_Arikan_GF4,1))-1
                g0 = kron(g0,g_Arikan_GF4);
                g0(g0==alpha_x_alpha) = alpha2;
                g0(g0==alpha_x_alpha2) = one;
                g0(g0==alpha2_x_alpha2) = alpha;
            end
            g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,2),'Arikan_GF4');
            
            Constellation_Mapping_Array = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            if(isnumeric(Lable_Type))
                
                date_string = strrep(sprintf('%s',datestr(now)),':','_');
                
                plot_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_gArikan_GF4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                file_name = sprintf('%s_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_gArikan_GF4',date_string,N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type);
                legend_name = sprintf('BIMLPCM %d%s gArikan GF4',M,Constellation_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
                
            else
            
                plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gArikan_GF4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gArikan_GF4',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type);
                legend_name = sprintf('BIMLPCM %d%s %s gArikan GF4',M,Constellation_Type,Lable_Type);
                title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
            
            end
            
        end
        
    end  
    
    SNR_Vec_dB_EsN0 = 10.*log10(SNR_Vec);
        
    SNR_Vec_dB_EbN0 = SNR_Vec_dB_EsN0 - 10*log10(m) - 10*log10(R);
    
    BLER_Axis_Buffer = 10^-1;
    Min_BLER = BLER_Axis_Buffer*BLER_Max;
    
    for i=1:1:size(SNR_Vec,2)

        if(sum(is_high_SNR)==size(is_high_SNR,2))
            break
        end
        
        if(R<1)
                            
            [Estimated_Bhattacharyya(i,:),SER_Vec(i,:),Is_Frozen_Bit_Index_Vec(i,:),~,Genie_Aided_LLR_Avg_Vec(i,:),Genie_Aided_LLR_Var_Vec(i,:),Estimated_SER_Vec(i,:),~,~,~,~,~] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,true,g,g0,Decode_Type,1,0,Constellation_Type,m,m_BIPCM,m_SCL,Constellation_Mapping_Array,N,R,BLER_Max,SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_Vec(i,:));
                        
            Plot_Index_SER(m,Estimated_Bhattacharyya(i,:),SER_Vec(i,:),Estimated_SER_Vec(i,:),Genie_Aided_LLR_Avg_Vec(i,:),sqrt(Genie_Aided_LLR_Var_Vec(i,:)),Is_Frozen_Bit_Index_Vec(i,:),title_name,legend_name,plot_path,sprintf('%s SNR=%s',plot_name,strrep(sprintf('%.3f',SNR_Vec_dB_EbN0(1,i)),'.','_')),is_visible);
            
        else
            
            Is_Frozen_Bit_Index_Vec(i,:) = ones(1,N);
            
        end

        Is_Frozen_Bit_Index_Vec(i,:) = Is_Frozen_Bit_Index_Vec(i,:) == ones(1,N);
        
        for j=1:1:size(L_List,2)
            
            if(is_high_SNR(1,j))
                continue
            end
            
            if(L_List(j)==1)
                
                Is_Frozen_Bit_Index_List_Vec(j,i,:) = Is_Frozen_Bit_Index_Vec(i,:);
                
                Is_Frozen_Bit_Index_temp_Vec(1:N) = Is_Frozen_Bit_Index_List_Vec(j,i,:);
                
                [~,~,~,~,~,~,~,BLER(1,i),BER(1,i),First_Symbol(1,i,:),LLR_Avg(1,i,:),LLR_Var(1,i,:)] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,false,g,g0,Decode_Type,L_List(j),GA_CRC_Length_List(j),Constellation_Type,m,m_BIPCM,m_SCL,Constellation_Mapping_Array,N,R,BLER_Max,SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_temp_Vec);
                
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
                
                [~,~,~,~,~,~,~,BLER(j,i),BER(j,i),First_Symbol(j,i,:),LLR_Avg(j,i,:),LLR_Var(j,i,:)] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,false,g,g0,Decode_Type,L_List(j),GA_CRC_Length_List(j),Constellation_Type,m,m_BIPCM,m_SCL,Constellation_Mapping_Array,N,R,BLER_Max,SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_temp_Vec);
        
            end
            
            temp_First_Symbol(1,:) = First_Symbol(j,i,:);
            temp_LLR_Avg(1,:) = LLR_Avg(j,i,:);
            temp_LLR_Var(1,:) = LLR_Var(j,i,:);
            
            if(L_List(j)~=1)
                if(GA_CRC_Length_List(j)>=0)
                    if(or(m_SCL~=m,m==1))
                        plot_name_temp = sprintf('L=%d_GA_CRC_%d_%s',L_List(j),GA_CRC_Length_List(j),plot_name);
                        legend_name_temp = sprintf('L=%d(%d/%d) with GA CRC %d %s',L_List(j),m_SCL,m,GA_CRC_Length_List(j),legend_name);
                    else
                        plot_name_temp = sprintf('L=%d(%dof%d)_GA_CRC_%d_%s',L_List(j),m_SCL,m,GA_CRC_Length_List(j),plot_name);
                        legend_name_temp = sprintf('L=%d with GA CRC %d %s',L_List(j),GA_CRC_Length_List(j),legend_name);
                    end
                else
                    if(or(m_SCL~=m,m==1))
                        plot_name_temp = sprintf('L=%d_%s',L_List(j),plot_name);
                        legend_name_temp = sprintf('L=%d(%d/%d) %s',L_List(j),m_SCL,m,legend_name);
                    else
                        plot_name_temp = sprintf('L=%d(%dof%d)_%s',L_List(j),m_SCL,m,plot_name);
                        legend_name_temp = sprintf('L=%d %s',L_List(j),legend_name);
                    end
                end
            else
                plot_name_temp = plot_name;
                legend_name_temp = legend_name;
            end
            
            Plot_Simulation_Results(m,SNR_Type,SNR_Vec_dB_EbN0(1,1:i),BLER(j,1:i),Min_BLER,BER(j,1:i),temp_First_Symbol,temp_LLR_Avg,sqrt(temp_LLR_Var),Is_Frozen_Bit_Index_Vec(i,:),title_name,legend_name_temp,plot_path,sprintf('%s SNR=%s',plot_name_temp,strrep(sprintf('%.3f',SNR_Vec_dB_EbN0(1,i)),'.','_')),is_visible);

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
    
end