function [BLER,BER,Optimal_BLER,Optimal_BER] = Find_Optimal_Frozen_Bits(NumOfWorkers,SNR_Type,Decode_Type,L_List,L_Optimize_List,GA_CRC_Length_List,Constellation_Type,m,m_BIPCM,Lable_Type,N,R,Optimal_Frozen_Bits_Search_Percentage,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible)

    g_Arikan = [1,0;1,1];

    M = 2^m;
    
    SER_Vec = zeros(length(SNR_Vec),N);
    Is_Frozen_Bit_Index_Vec = zeros(length(SNR_Vec),N);
    Is_Frozen_Bit_Index_Vec = Is_Frozen_Bit_Index_Vec == ones(length(SNR_Vec),N);
    Is_Frozen_Bit_Index_List_Vec = zeros(length(L_List),length(SNR_Vec),N);
    Is_Frozen_Bit_Index_List_Diff_Vec = zeros(length(L_List),length(SNR_Vec));
    Genie_Aided_LLR_Avg_Vec = zeros(length(SNR_Vec),N);
    Genie_Aided_LLR_Var_Vec = zeros(length(SNR_Vec),N);
    Estimated_SER_Vec = zeros(length(SNR_Vec),N);
    BLER = zeros(length(L_List),length(SNR_Vec));
    BER = zeros(length(L_List),length(SNR_Vec));
    First_Bit = zeros(length(L_List),length(SNR_Vec),N);
    temp_First_Bit = zeros(1,N);
    LLR_Avg = zeros(length(SNR_Vec),N);
    LLR_Var = zeros(length(SNR_Vec),N);
    Is_Change_Frozen_Bits_Index_Vec = zeros(length(SNR_Vec)-1,N);
    Joint_SER_Vec = zeros(length(SNR_Vec),N);
    Optimal_BLER = zeros(length(L_List),length(SNR_Vec));
    Optimal_BER = zeros(length(L_List),length(SNR_Vec));
    Optimal_First_Bit = zeros(length(L_List),length(SNR_Vec),N);
    Optimal_LLR_Avg = zeros(length(SNR_Vec),N);
    Optimal_LLR_Var = zeros(length(SNR_Vec),N);
    Is_Optimal_Frozen_Bit_Index_Vec = zeros(length(SNR_Vec),N);
    Is_Optimal_Frozen_Bit_Index_Vec = Is_Optimal_Frozen_Bit_Index_Vec == ones(length(SNR_Vec),N);
    Is_Optimal_Frozen_Bit_Index_Diff_Vec = zeros(length(SNR_Vec),N);
    Estimated_Bhattacharyya = zeros(length(SNR_Vec),N);

    if(R>1)
        BLER = ones(length(L_List),length(SNR_Vec));
        BER = ones(length(L_List),length(SNR_Vec));
        return;
    end

    if(M==2)
                
        g = g_Arikan;
        for i=1:1:log2(N)-1
            g = kron(g,g_Arikan);
        end
        g = bitrevorder(g);
        
        g0 = zeros(0,0);
        
        Constellation_Mapping_Array = [-1,1];
        
        plot_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BPSK_gArikan',N,strrep(sprintf('%.3f',R),'.','_'));
        file_name = sprintf('N=%d[bits]_R=%s[bits_per_channel_use]_BPSK_gArikan',N,strrep(sprintf('%.3f',R),'.','_'));
        legend_name = sprintf('BPSK gArikan');
        title_name = sprintf('N=%d[bits] R=%s[Information bits/Total bits]',N,sprintf('%.3f',R));
        
    elseif(mod(log(M)/log(4),1)==0)
                
        if(strcmp(Decode_Type,'MLPC'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = zeros(0,0);
            
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
            
            g0 = zeros(0,0);
            
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
            
            %g0 = zeros(0,0);
            
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
            
        elseif(strcmp(Decode_Type,'Compound_BIPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = g_Arikan;
            g0 = bitrevorder(g0);
            g0 = kron(eye(m/2),g0);

            %g0 = g_Arikan;
            %for i=1:1:log2(m)-1
            %    g0 = kron(g0,g_Arikan);
            %end
            %g0 = bitrevorder(g0);
            
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
            
        end
        
    end  
    
    SNR_Vec_dB_EsN0 = 10.*log10(SNR_Vec);
        
    SNR_Vec_dB_EbN0 = SNR_Vec_dB_EsN0 - 10*log10(m) - 10*log10(R);
    
    BLER_Axis_Buffer = 10^-1;
    Min_BLER = BLER_Axis_Buffer*BLER_Max;
    
    for i=1:1:length(SNR_Vec)

        if(R<1)
                            
            [Estimated_Bhattacharyya(i,:),SER_Vec(i,:),Is_Frozen_Bit_Index_Vec(i,:),~,Genie_Aided_LLR_Avg_Vec(i,:),Genie_Aided_LLR_Var_Vec(i,:),Estimated_SER_Vec(i,:),~,~,~,~,~] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,true,g,g0,Decode_Type,1,0,Constellation_Type,m,m_BIPCM,0,Constellation_Mapping_Array,N,R,BLER_Max,SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_Vec(i,:));
                        
            %Plot_Index_SER(m,Estimated_Bhattacharyya(i,:),SER_Vec(i,:),Estimated_SER_Vec(i,:),Genie_Aided_LLR_Avg_Vec(i,:),sqrt(Genie_Aided_LLR_Var_Vec(i,:)),Is_Frozen_Bit_Index_Vec(i,:),title_name,legend_name,plot_path,sprintf('%s SNR=%s',plot_name,strrep(sprintf('%.3f',SNR_Vec_dB_EbN0(1,i)),'.','_')),is_visible);
            
        else
            
            Is_Frozen_Bit_Index_Vec(i,:) = ones(1,N);
            
        end

        Is_Frozen_Bit_Index_Vec(i,:) = Is_Frozen_Bit_Index_Vec(i,:) == ones(1,N);
        
        for j=1:1:length(L_List)
            
            if(L_List(j)~=1)
                if(GA_CRC_Length_List(j)>=0)
                    if(or(m_SCL~=m,m==1))
                        plot_name_temp = sprintf('L=%d_GA_CRC_%d_%s',L_List(j),GA_CRC_Length_List(j),plot_name);
                        legend_name_temp = sprintf('L=%d(%d/%d) with GA CRC %d %s',L,m_SCL,m,GA_CRC_Length_List(j),legend_name);
                    else
                        plot_name_temp = sprintf('L=%d(%dof%d)_GA_CRC_%d_%s',L_List(j),m_SCL,m,GA_CRC_Length_List(j),plot_name);
                        legend_name_temp = sprintf('L=%d with GA CRC %d %s',L,GA_CRC_Length_List(j),legend_name);
                    end
                else
                    if(or(m_SCL~=m,m==1))
                        plot_name_temp = sprintf('L=%d_%s',L_List(j),plot_name);
                        legend_name_temp = sprintf('L=%d(%d/%d) %s',L,m_SCL,m,legend_name);
                    else
                        plot_name_temp = sprintf('L=%d(%dof%d)_%s',L_List(j),m_SCL,m,plot_name);
                        legend_name_temp = sprintf('L=%d %s',L,legend_name);
                    end
                end
            else
                plot_name_temp = plot_name;
                legend_name_temp = legend_name;
            end
            
            if(L_List(j)==1)
                
                Is_Frozen_Bit_Index_List_Vec(j,i,:) = Is_Frozen_Bit_Index_Vec(i,:);
                
                [~,~,~,~,~,~,~,BLER(1,i),BER(1,i),First_Bit(1,i,:),LLR_Avg(i,:),LLR_Var(i,:)] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,false,g,g0,Decode_Type,L_List(j),GA_CRC_Length_List(j),Constellation_Type,m,m_BIPCM,0,Constellation_Mapping_Array,N,R,BLER_Max,SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_Vec(i,:));
                
                temp_First_Bit(1,:) = First_Bit(1,i,:);
                
                Plot_Simulation_Results(m,SNR_Type,SNR_Vec_dB_EbN0(1,1:i),BLER(1,1:i),Min_BLER,BER(1,1:i),temp_First_Bit,LLR_Avg(i,:),sqrt(LLR_Var(i,:)),Is_Frozen_Bit_Index_Vec(i,:),title_name,legend_name_temp,plot_path,sprintf('%s SNR=%s',plot_name_temp,strrep(sprintf('%.3f',SNR_Vec_dB_EbN0(1,i)),'.','_')),is_visible);
                
                if(i>1)
                    Is_Change_Frozen_Bits_Index_Vec(i-1,:) = Is_Frozen_Bit_Index_Vec(i-1,:)~=Is_Frozen_Bit_Index_Vec(i,:);
                    sprintf('%d frozen bits has changed',sum(Is_Change_Frozen_Bits_Index_Vec(i-1,:))/2)
                end
                
                [Sorted_SER_Vec,Sorted_SER_Index_Vec] = sort(wrev(SER_Vec(i,:)));
                Sorted_SER_Index_Vec = length(SER_Vec) - Sorted_SER_Index_Vec + 1;

                temp_R = floor((R+Optimal_Frozen_Bits_Search_Percentage)*N)/N;
                
                Is_Optimal_Frozen_Bit_Index_Vec(i,:) = zeros(1,N);
                for Index=Sorted_SER_Index_Vec(1,temp_R*N+1+(GA_CRC_Length>0)*GA_CRC_Length:N)
                    Is_Optimal_Frozen_Bit_Index_Vec(i,Index) = 1;
                end
                Is_Optimal_Frozen_Bit_Index_Vec(i,:) = Is_Optimal_Frozen_Bit_Index_Vec(i,:) == ones(1,N);
                
                [~,~,~,~,~,~,~,Optimal_BLER(1,i),Optimal_BER(1,i),First_Bit(1,i,:),LLR_Avg(i,:),LLR_Var(i,:)] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,false,g,g0,Decode_Type,L_List(j),GA_CRC_Length_List(j),Constellation_Type,m,m_BIPCM,0,Constellation_Mapping_Array,N,temp_R,BLER_Max,SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Optimal_Frozen_Bit_Index_Vec(i,:));
                
                Joint_SER_Vec(i,First_Bit(1,i,:)~=0) = First_Bit(1,i,First_Bit(1,i,:)~=0);
                Joint_SER_Vec(i,First_Bit(1,i,:)==0) = SER_Vec(i,First_Bit(1,i,:)==0);
                
                [Sorted_SER_Vec,Sorted_SER_Index_Vec] = sort(wrev(Joint_SER_Vec(i,:)));
                Sorted_SER_Index_Vec = length(Joint_SER_Vec) - Sorted_SER_Index_Vec + 1;

                Is_Optimal_Frozen_Bit_Index_Vec(i,:) = zeros(1,N);
                for Index=Sorted_SER_Index_Vec(1,R*N+1+(GA_CRC_Length>0)*GA_CRC_Length:N)
                    Is_Optimal_Frozen_Bit_Index_Vec(i,Index) = 1;
                end
                Is_Optimal_Frozen_Bit_Index_Vec(i,:) = Is_Optimal_Frozen_Bit_Index_Vec(i,:) == ones(1,N);
                
                Is_Optimal_Frozen_Bit_Index_Diff_Vec(i,:) = Is_Optimal_Frozen_Bit_Index_Vec(i,:)~=Is_Frozen_Bit_Index_Vec(i,:);
                
                sprintf('optimal frozen bits: %d frozen bits has changed',sum(Is_Optimal_Frozen_Bit_Index_Diff_Vec(i,:))/2)
                
                [~,~,~,~,~,~,~,Optimal_BLER(1,i),Optimal_BER(1,i),Optimal_First_Bit(1,i,:),Optimal_LLR_Avg(i,:),Optimal_LLR_Var(i,:)] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,false,g,g0,Decode_Type,L_List(j),GA_CRC_Length_List(j),Constellation_Type,m,m_BIPCM,0,Constellation_Mapping_Array,N,R,BLER_Max,SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Optimal_Frozen_Bit_Index_Vec(i,:));
                
                temp_First_Bit(1,:) = Optimal_First_Bit(1,i,:);
                
                Plot_Simulation_Results(m,1,1,0,SNR_Type,SNR_Vec_dB_EbN0(1,1:i),Optimal_BLER(1,1:i),Min_BLER,Optimal_BER(1,1:i),temp_First_Bit,LLR_Avg(i,:),sqrt(LLR_Var(i,:)),Is_Frozen_Bit_Index_Vec(i,:),title_name,legend_name,plot_path,sprintf('Optimal Frozen Bits %s SNR=%s',plot_name,strrep(sprintf('%.3f',SNR_Vec_dB_EbN0(1,i)),'.','_')),is_visible);
                
            else
             
                [Sorted_SER_Vec,Sorted_SER_Index_Vec] = sort(wrev(SER_Vec(i,:)));
                Sorted_SER_Index_Vec = length(SER_Vec(i,:)) - Sorted_SER_Index_Vec + 1;

                [Is_Frozen_Bit_Index_List_Vec(j,i,:),Is_Frozen_Bit_Index_List_Diff_Vec(j,i)] = Find_Frozen_Bits_List(Sorted_SER_Vec,Sorted_SER_Index_Vec,R,L_Optimize_List(j));
            
                [~,~,~,~,~,~,~,BLER(j,i),BER(j,i),First_Bit(j,i,:),~,~] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,false,g,g0,Decode_Type,L_List(j),GA_CRC_Length_List(j),Constellation_Type,m,m_BIPCM,0,Constellation_Mapping_Array,N,R,BLER_Max,SNR_Vec(i),Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_List_Vec(j,i,:));
        
            end
            
        end
        
        try
            save(fullfile(file_path,file_name),'SNR_Vec_dB_EbN0','SNR_Vec_dB_EsN0','SER_Vec','Genie_Aided_LLR_Avg_Vec','Genie_Aided_LLR_Var_Vec','Estimated_SER_Vec','BLER','BER','Optimal_BLER','Optimal_BER','First_Bit','LLR_Avg','LLR_Var','Optimal_First_Bit','Optimal_LLR_Avg','Optimal_LLR_Var','Is_Frozen_Bit_Index_Vec','Is_Change_Frozen_Bits_Index_Vec','Is_Frozen_Bit_Index_List_Vec','Is_Frozen_Bit_Index_List_Diff_Vec','Is_Optimal_Frozen_Bit_Index_Vec');
        end
        
        if(and(BLER(1,i)<BLER_Max*10,i<length(SNR_Vec)))
            for j=i+1:1:length(SNR_Vec)
                Is_Frozen_Bit_Index_Vec(j,:) = Is_Frozen_Bit_Index_Vec(i,:);
                BLER(:,j) = zeros(size(BLER,1),1);
            end
            break;
        end
        
    end
    
    sprintf('%d frozen bits has changed\n',sum(Is_Change_Frozen_Bits_Index_Vec,2)./2)
    
end