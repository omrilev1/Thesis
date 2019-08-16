function [] = N_m_PDF(NumOfWorkers,Decode_Type,Constellation_Type,m,m_BIPCM,Constellation_Mapping_Arrays,Lable_Type_Name,R,SNR_Vec,BLER_Max,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,legend_name,is_visible)

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
    
    N = m;
    
    M = 2^m;
    
    if(strcmp(Decode_Type,'gRS4'))
        U_N = N/2;
        m_GF = m/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'MLPC_gRS4'))
        U_N = N/2;
        m_GF = m/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIPCM_gRS4'))
        U_N = N/2;
        m_GF = m/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIMLPCM_gRS4'))
        U_N = N/2;
        m_GF = m/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'gArikan_GF4'))
        U_N = N/2;
        m_GF = m/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
        U_N = N/2;
        m_GF = m/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIPCM_gArikan_GF4'))
        U_N = N/2;
        m_GF = m/2;
        Is_GF4 = true;
    elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4'))
        U_N = N/2;
        m_GF = m/2;
        Is_GF4 = true;
    else %gArikan
        U_N = N;
        m_GF = m;
        Is_GF4 = false;
    end
    
    if(M==2)
                
        return

    elseif(mod(log(M)/log(4),1)==0)
        
        if(strcmp(Decode_Type,'MLPC'))
            
            g = 1;
            g0 = 1;
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_MLPC_%d%s_%s_gArikan',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_MLPC_%d%s_gArikan',date_string,N,M,Constellation_Type);
            title_name = sprintf('m=%d[bits]',m);
            
        elseif(strcmp(Decode_Type,'MLPC_gRS4'))
            
            g = 1;
            g0 = 1;
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_MLPC_%d%s_%s_gRS4',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_MLPC_%d%s_gRS4',date_string,N,M,Constellation_Type);
            
        elseif(strcmp(Decode_Type,'MLPC_gArikan_GF4'))
            
            g = 1;
            g0 = 1;
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_MLPC_%d%s_%s_gArikan_GF4',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_MLPC_%d%s_gArikan_GF4',date_string,N,M,Constellation_Type);
            
        elseif(strcmp(Decode_Type,'Compound_MLPC'))
            
%             g = g_Arikan;
%             for i=1:1:log2(N/m)-1
%                 g = kron(g,g_Arikan);
%             end
%             g = bitrevorder(g);
            
            g = 1;

            g0 = g_Arikan;
            for i=1:1:log2(m)-1
                g0 = kron(g0,g_Arikan);
            end
            g0 = bitrevorder(g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_Compound_MLPC_%d%s_%s_gArikan',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_Compound_MLPC_%d%s_gArikan',date_string,N,M,Constellation_Type);
            title_name = sprintf('m=%d[bits]',m);
            
        elseif(strcmp(Decode_Type,'Separated_BIPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = zeros(0,0);
             
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_Separated_BIPCM_%d%s_%s_gArikan',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_Separated_BIPCM_%d%s_gArikan',date_string,N,M,Constellation_Type);
            title_name = sprintf('m=%d[bits]',m);
            
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
                        
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_BIPCM_%d%s_%s_gArikan',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_BIPCM_%d%s_gArikan',date_string,N,M,Constellation_Type);
            title_name = sprintf('m=%d[bits]',m);
        
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
                   
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_BIPCM_%d%s_%s_gRS4',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_BIPCM_%d%s_gRS4',date_string,N,M,Constellation_Type);
            title_name = sprintf('m=%d[bits]',m);    
         
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
                     
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_BIPCM_%d%s_%s_gArikan_GF4',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_BIPCM_%d%s_gArikan_GF4',date_string,N,M,Constellation_Type);
            title_name = sprintf('m=%d[bits]',m);    
            
        elseif(strcmp(Decode_Type,'Compound_BIPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = g_Arikan;
            g0 = bitrevorder(g0);
            g0 = kron(eye(m/2),g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_Compound_BIPCM_%d%s_%s_gArikan',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_Compound_BIPCM_%d%s_gArikan',date_string,N,M,Constellation_Type);
            title_name = sprintf('m=%d[bits]',m);
            
         elseif(strcmp(Decode_Type,'BIMLPCM'))
            
%             g = g_Arikan;
%             for i=1:1:log2(N/m)-1
%                 g = kron(g,g_Arikan);
%             end
%             g = bitrevorder(g);

            g = 1;
            
            g0 = g_Arikan;
            for i=1:1:log2(m)-1
                g0 = kron(g0,g_Arikan);
            end
            g0 = bitrevorder(g0);
                        
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_BIMLPCM_%d%s_%s_gArikan',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_BIMLPCM_%d%s_gArikan',date_string,N,M,Constellation_Type);
            title_name = sprintf('m=%d[bits]',m);
           
        elseif(strcmp(Decode_Type,'BIMLPCM_gRS4'))
            
%             g = g_RS4;
%             for i=1:1:log2(N/m)/log2(size(g_RS4,1))-1
%                 g = kron(g,g_RS4);
%                 g(g==alpha_x_alpha) = alpha2;
%                 g(g==alpha_x_alpha2) = one;
%                 g(g==alpha2_x_alpha2) = alpha;
%             end
%             g = Convert_Matrix_GF4_to_Matrix(digitrevorder(g,4),'RS4');

            g = 1;

            g0 = g_RS4;
            for i=1:1:log2(m)/log2(2*size(g_RS4,1))-1
                g0 = kron(g0,g_RS4);
                g0(g0==alpha_x_alpha) = alpha2;
                g0(g0==alpha_x_alpha2) = one;
                g0(g0==alpha2_x_alpha2) = alpha;
            end
            g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,4),'RS4');
                   
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_BIMLPCM_%d%s_%s_gRS4',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_BIMLPCM_%d%s_gRS4',date_string,N,M,Constellation_Type);
            title_name = sprintf('m=%d[bits]',m);
            
        elseif(strcmp(Decode_Type,'BIMLPCM_gArikan_GF4'))
            
%             g = g_Arikan_GF4;
%             for i=1:1:log2(N/m)/log2(size(g_Arikan_GF4,1))-1
%                 g = kron(g,g_Arikan_GF4);
%                 g(g==alpha_x_alpha) = alpha2;
%                 g(g==alpha_x_alpha2) = one;
%                 g(g==alpha2_x_alpha2) = alpha;
%             end
%             g = Convert_Matrix_GF4_to_Matrix(digitrevorder(g,2),'Arikan_GF4');

            g = 1;

            g0 = g_Arikan_GF4;
            for i=1:1:log2(m)/log2(2*size(g_Arikan_GF4,1))-1
                g0 = kron(g0,g_Arikan_GF4);
                g0(g0==alpha_x_alpha) = alpha2;
                g0(g0==alpha_x_alpha2) = one;
                g0(g0==alpha2_x_alpha2) = alpha;
            end
            g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,2),'Arikan_GF4');
                  
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            file_name = sprintf('%s_Same_Type_Labels_N=%d[bits]_BIMLPCM_%d%s_%s_gArikan_GF4',date_string,N,M,Constellation_Type,Lable_Type_Name);
            plot_name = sprintf('%s_N=%d[bits]_BIMLPCM_%d%s_gArikan_GF4',date_string,N,M,Constellation_Type);
            title_name = sprintf('m=%d[bits]',m);
            
        end
        
    end  
    
    SNR_Vec_dB_EsN0 = 10.*log10(SNR_Vec);
    SNR_Vec_dB_EbN0 = SNR_Vec_dB_EsN0 - 10*log10(m) - 10*log10(R);
    
    Num_Of_Iterations = 1/BLER_Max;
    Num_Of_Bins = floor(Num_Of_Iterations/100);
    Num_Of_STD = 10;
    LLR_PDF_Vec = zeros(size(SNR_Vec,2),size(Constellation_Mapping_Arrays,1),m_GF,Num_Of_Bins);
    temp_LLR_PDF_Vec = zeros(size(Constellation_Mapping_Arrays,1),m_GF,Num_Of_Bins);
    LLR_Bin_Border_Vec = zeros(size(SNR_Vec,2),size(Constellation_Mapping_Arrays,1),m_GF,Num_Of_Bins);
    temp_LLR_Bin_Border_Vec = zeros(size(Constellation_Mapping_Arrays,1),m_GF,Num_Of_Bins);
    LLR_Bin_Border_Vec_temp =zeros(m_GF,Num_Of_Bins);
    
    try
        save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','Sum_SER_Vec');
    end
    
    for i=1:1:size(SNR_Vec,2)
        
        for j=1:1:size(Constellation_Mapping_Arrays,1)

            Constellation_Mapping_Array = Constellation_Mapping_Arrays(j,:);

            [~,~,~,~,Genie_Aided_LLR_Avg_Vec,Genie_Aided_LLR_Var_Vec,~,~,~,~,~,~] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,true,g,g0,Decode_Type,1,0,Constellation_Type,m,m_BIPCM,0,Constellation_Mapping_Array,N,(N-1)/N,BLER_Max,SNR_Vec(1,i),Min_Num_Of_Errors_Sufficient_Statistics,zeros(1,N));
            
            for k=1:1:m_GF
                LLR_Bin_Border_Vec(i,j,k,:) = Genie_Aided_LLR_Avg_Vec(1,k)-Num_Of_STD*sqrt(Genie_Aided_LLR_Var_Vec(1,k)):2*Num_Of_STD*sqrt(Genie_Aided_LLR_Var_Vec(1,k))/(Num_Of_Bins-1):Genie_Aided_LLR_Avg_Vec(1,k)+Num_Of_STD*sqrt(Genie_Aided_LLR_Var_Vec(1,k));
            end
            
            LLR_Bin_Border_Vec_temp(:,:) = LLR_Bin_Border_Vec(i,j,:,:);
            
            [LLR_PDF_Vec(i,j,:,:)] = Simulate_Constellation_Mapping_PDF(NumOfWorkers,g,g0,Decode_Type,Constellation_Type,m,m_BIPCM,Constellation_Mapping_Array,Num_Of_Iterations,SNR_Vec(1,i),LLR_Bin_Border_Vec_temp);

            try
                save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','LLR_PDF_Vec');
            end

        end
        
        temp_LLR_PDF_Vec(:,:,:) = LLR_PDF_Vec(i,:,:,:);
        temp_LLR_Bin_Border_Vec(:,:,:) = LLR_Bin_Border_Vec(i,:,:,:);
        
        Plot_LLR_PDF(m_GF,temp_LLR_Bin_Border_Vec,temp_LLR_PDF_Vec,title_name,legend_name,plot_path,sprintf('%s_%s',plot_name,strrep(sprintf('%.3f',SNR_Vec_dB_EbN0(1,i)),'.','_')),is_visible);
        
    end
    
    
end