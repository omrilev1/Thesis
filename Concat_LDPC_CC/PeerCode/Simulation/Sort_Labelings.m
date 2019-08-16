function [Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound,Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound,Sorted_By_Union_Bound] = Sort_Labelings(NumOfWorkers,Decode_Type,Constellation_Type,m,m_BIPCM,Lable_Type,N,R,SNR_Vec,Estimted_BLER,Min_Num_Of_Errors_Sufficient_Statistics,file_path,plot_path,is_visible,is_sorted_by_N_m,is_sorted_by_Union_Bound)

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
    
    M = 2^m;
    
    Num_Of_Iterations_Print_State = 100;
    
%     SNR_Constellation_Polar_Code = SNR_Vec(1,1);
    SNR_Constellation_Polar_Code = SNR_Vec(1,end);
    SNR_Constellation_Union_Bound = SNR_Vec(1,end);
    
%     Estimted_BLER = Estimted_BLER_Vec(1,1);
    Estimted_BLER_Sorted_By_Union_Bound = 10^0;
%     Estimted_BLER_High_SNR = Estimted_BLER_Vec(1,end);

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
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = 1;
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gArikan_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);

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
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gRS4_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);

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
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_MLPC_%d%s_%s_gArikan_GF4_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);
            
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
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);

            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_Compound_MLPC_%d%s_%s_gArikan_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);
            
        elseif(strcmp(Decode_Type,'Separated_BIPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = 1;
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
 
            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_Separated_BIPCM_%d%s_%s_gArikan_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);
            
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
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_%s_gArikan_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);
            
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
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_%s_gRS4_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);
            
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
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_BIPCM_%d%s_%s_gArikan_GF4_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);
            
        elseif(strcmp(Decode_Type,'Compound_BIPCM'))
            
            g = g_Arikan;
            for i=1:1:log2(N/m)-1
                g = kron(g,g_Arikan);
            end
            g = bitrevorder(g);
            
            g0 = g_Arikan;
            g0 = bitrevorder(g0);
            g0 = kron(eye(m/2),g0);
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);

            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_Compound_BIPCM_%d%s_%s_gArikan_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);
            
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
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gArikan_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);
            
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
            for i=1:1:log2(m)/log2(2*size(g_RS4,1))-1
                g0 = kron(g0,g_RS4);
                g0(g0==alpha_x_alpha) = alpha2;
                g0(g0==alpha_x_alpha2) = one;
                g0(g0==alpha2_x_alpha2) = alpha;
            end
            g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,4),'RS4');
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gRS4_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);
            
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
            for i=1:1:log2(m)/log2(2*size(g_Arikan_GF4,1))-1
                g0 = kron(g0,g_Arikan_GF4);
                g0(g0==alpha_x_alpha) = alpha2;
                g0(g0==alpha_x_alpha2) = one;
                g0(g0==alpha2_x_alpha2) = alpha;
            end
            g0 = Convert_Matrix_GF4_to_Matrix(digitrevorder(g0,2),'Arikan_GF4');
            
            Constellation_Mapping_Arrays = Make_Constellation(Constellation_Type,m,Lable_Type,g0);
            
            [~,temp_Lable_Type] = strtok(Lable_Type,'-');
            if(isempty(temp_Lable_Type))
                temp_Lable_Type = Lable_Type;
            else
                temp_Lable_Type = temp_Lable_Type(2:end);
            end
            Regular_Labeling = Make_Constellation(Constellation_Type,m,temp_Lable_Type,g0);
            
            date_string = strrep(sprintf('%s',datestr(now)),':','_');
            
            file_name = sprintf('Same_Type_Labels_N=%d[bits]_R=%s[bits_per_channel_use]_BIMLPCM_%d%s_%s_gArikan_GF4_%s',N,strrep(sprintf('%.3f',R),'.','_'),M,Constellation_Type,Lable_Type,date_string);
            
        end
        
    end  
    
    Regular_Labeling_Index = find(~sum(abs(bsxfun(@minus,Constellation_Mapping_Arrays,Regular_Labeling)),2),1);
    if(isempty(Regular_Labeling_Index))
        Regular_Labeling_Index = size(Constellation_Mapping_Arrays,1)+1;
        Constellation_Mapping_Arrays(Regular_Labeling_Index,:) = Regular_Labeling;
    end
    
    SER_Vec = zeros(size(Constellation_Mapping_Arrays,1),U_N);
    Is_Frozen_Bit_Index_Vec = zeros(size(Constellation_Mapping_Arrays,1),N);
    Estimated_Is_Frozen_Bit_Index_Vec = zeros(size(Constellation_Mapping_Arrays,1),N);
    Genie_Aided_LLR_Avg_Vec = zeros(size(Constellation_Mapping_Arrays,1),U_N);
    Genie_Aided_LLR_Var_Vec = zeros(size(Constellation_Mapping_Arrays,1),U_N);
    Estimated_SER_Vec = zeros(size(Constellation_Mapping_Arrays,1),U_N);
    Estimated_Bhattacharyya = zeros(size(Constellation_Mapping_Arrays,1),U_N);
    
    Num_Of_Iterations = floor((Min_Num_Of_Errors_Sufficient_Statistics*(1/Estimted_BLER)));
    %Num_Of_Bins = floor(Num_Of_Iterations/100);
    %Num_Of_Bins = floor(Num_Of_Iterations/10);
    Num_Of_Bins = 2;
%     Num_Of_STD = 10;
    
%     N_m_PDF_Genie_Aided_LLR_Avg_Vec = zeros(size(Constellation_Mapping_Arrays,1),m);
%     N_m_PDF_Genie_Aided_LLR_Var_Vec = zeros(size(Constellation_Mapping_Arrays,1),m);

    LLR_PDF_Vec = zeros(size(Constellation_Mapping_Arrays,1),m_GF,Num_Of_Bins);
    LLR_Bin_Border_Vec = zeros(size(Constellation_Mapping_Arrays,1),m_GF,Num_Of_Bins);
%         LLR_Bin_Border_Vec_temp =zeros(m_GF,Num_Of_Bins);
    Estimated_BER = zeros(size(Constellation_Mapping_Arrays,1),m_GF);

%     temp_LLR_PDF_Vec = zeros(1,Num_Of_Bins);
%     temp_LLR_Bin_Border_Vec = zeros(1,Num_Of_Bins);
    
    Sorted_By_Union_Bound = zeros(0,0);
    Sorted_By_Estimated_Constellation_Polar_Code_LLR = zeros(0,0);
    
    d_min_square = zeros(size(Constellation_Mapping_Arrays,1),m_GF);
    Es = zeros(size(Constellation_Mapping_Arrays,1),m_GF);
    d_min_square_to_Es = zeros(size(Constellation_Mapping_Arrays,1),m_GF);
    
    try
        save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays');
    end
       
    for i=1:1:size(Constellation_Mapping_Arrays,1)
        
        if(is_sorted_by_Union_Bound)
            
            [Estimated_Bhattacharyya(i,:),SER_Vec(i,:),Is_Frozen_Bit_Index_Vec(i,:),Estimated_Is_Frozen_Bit_Index_Vec(i,:),Genie_Aided_LLR_Avg_Vec(i,:),Genie_Aided_LLR_Var_Vec(i,:),Estimated_SER_Vec(i,:),~,~,~,~,~] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,true,g,g0,Decode_Type,1,0,Constellation_Type,m,m_BIPCM,0,Constellation_Mapping_Arrays(i,:),N,R,Estimted_BLER_Sorted_By_Union_Bound,SNR_Constellation_Union_Bound,Min_Num_Of_Errors_Sufficient_Statistics,Is_Frozen_Bit_Index_Vec(i,:));
        
        end
        
        if(is_sorted_by_N_m)
            
%             [~,~,~,~,N_m_PDF_Genie_Aided_LLR_Avg_Vec(i,:),N_m_PDF_Genie_Aided_LLR_Var_Vec(i,:),~,~,~,~,~,~] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,true,g0,g0,Decode_Type,1,0,Constellation_Type,m,m_BIPCM,0,Constellation_Mapping_Arrays(i,:),m,(m-1)/m,Estimted_BLER,SNR_Constellation_Polar_Code,Min_Num_Of_Errors_Sufficient_Statistics,zeros(1,N));
% 
%             for j=1:1:m
%                 LLR_Bin_Border_Vec(i,j,:) = N_m_PDF_Genie_Aided_LLR_Avg_Vec(i,j)-Num_Of_STD*sqrt(N_m_PDF_Genie_Aided_LLR_Var_Vec(i,j)):2*Num_Of_STD*sqrt(N_m_PDF_Genie_Aided_LLR_Var_Vec(i,j))/(Num_Of_Bins-1):N_m_PDF_Genie_Aided_LLR_Avg_Vec(i,j)+Num_Of_STD*sqrt(N_m_PDF_Genie_Aided_LLR_Var_Vec(i,j));
%             end
% 
%             LLR_Bin_Border_Vec_temp(:,:) = LLR_Bin_Border_Vec(i,:,:);
% 
%             [LLR_PDF_Vec(i,:,:)] = Simulate_Constellation_Mapping_PDF(NumOfWorkers,g0,g0,Decode_Type,Constellation_Type,m,m_BIPCM,Constellation_Mapping_Arrays(i,:),Num_Of_Iterations,Low_SNR,LLR_Bin_Border_Vec_temp);
        
            [LLR_PDF_Vec(i,:,:)] = Simulate_Constellation_Mapping_PDF(NumOfWorkers,g0,g0,Decode_Type,Constellation_Type,m,m_BIPCM,Constellation_Mapping_Arrays(i,:),Num_Of_Iterations,SNR_Constellation_Polar_Code,[-1,1]);
%             Estimated_BER(i,:) = LLR_PDF_Vec(i,:,1);

        end
        %Plot_Index_SER(m,Estimated_Bhattacharyya(i,:),SER_Vec(i,:),Estimated_SER_Vec(i,:),Genie_Aided_LLR_Avg_Vec(i,:),sqrt(Genie_Aided_LLR_Var_Vec(i,:)),Is_Frozen_Bit_Index_Vec(i,:),'SER','',plot_path,'plot',is_visible);
        %Plot_Index_SER(m,Estimated_Bhattacharyya(i,:),Estimated_SER_Vec(i,:),Estimated_SER_Vec(i,:),Genie_Aided_LLR_Avg_Vec(i,:),sqrt(Genie_Aided_LLR_Var_Vec(i,:)),Estimated_Is_Frozen_Bit_Index_Vec(i,:),'Estimated SER','',plot_path,'plot',is_visible);
        
        if(mod(i,Num_Of_Iterations_Print_State)==0)
            try
                if(and(is_sorted_by_N_m,is_sorted_by_Union_Bound))
                    save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','Estimated_SER_Vec','Estimated_Is_Frozen_Bit_Index_Vec','LLR_Bin_Border_Vec','LLR_PDF_Vec','Estimated_BER');
                elseif(and(is_sorted_by_N_m,~is_sorted_by_Union_Bound))
                    save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','LLR_Bin_Border_Vec','LLR_PDF_Vec','Estimated_BER');
                elseif(and(~is_sorted_by_N_m,is_sorted_by_Union_Bound))
                    save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','Estimated_SER_Vec','Estimated_Is_Frozen_Bit_Index_Vec');
                end
            end
            sprintf('finish label %d (out of %d labels)',i,size(Constellation_Mapping_Arrays,1))
        end
        
    end
    
    if(is_sorted_by_Union_Bound)
        
        if(Is_GF4)
            Union_Bound_Vec = sum(((~Estimated_Is_Frozen_Bit_Index_Vec(:,1:2:end)).*Estimated_SER_Vec).');
        else
            Union_Bound_Vec = sum(((~Estimated_Is_Frozen_Bit_Index_Vec).*Estimated_SER_Vec).');
        end
        [~,Sorted_Union_Bound_Index_Vec] = sort(Union_Bound_Vec);
        Sorted_By_Union_Bound = Constellation_Mapping_Arrays(Sorted_Union_Bound_Index_Vec,:);

    end
    
    if(is_sorted_by_N_m)
        
%         [~,~,Is_Frozen_Bit_Index_Vec_Regular_Labeling,~,~,~,~,~,~,~,~,~] = Find_Frozen_Bits_Or_Simulate(NumOfWorkers,true,g,g0,Decode_Type,1,0,Constellation_Type,m,m_BIPCM,0,Constellation_Mapping_Arrays(Regular_Labeling_Index,:),N,R,Estimted_BLER,SNR_Constellation_Polar_Code,Min_Num_Of_Errors_Sufficient_Statistics,zeros(1,N));
%         if(Is_GF4)
%             Estimated_Rates = (((N/m)-sum(reshape(Is_Frozen_Bit_Index_Vec_Regular_Labeling(1:2:end).',[],m_GF).',2))/(N/m)).';
%         else
%             Estimated_Rates = (((N/m)-sum(reshape(Is_Frozen_Bit_Index_Vec_Regular_Labeling.',[],m).',2))/(N/m)).';
%         end
%         for i=1:1:size(LLR_Bin_Border_Vec,1)
%             for j=1:1:m
%                 temp_LLR_Bin_Border_Vec(:) = LLR_Bin_Border_Vec(i,j,:);
%                 temp_LLR_PDF_Vec(:) = LLR_PDF_Vec(i,j,:);
%                 Estimated_BER(i,j) = sum(temp_LLR_PDF_Vec(temp_LLR_Bin_Border_Vec<0));
%             end
%         end
        for i=1:1:size(LLR_Bin_Border_Vec,1)
            for j=1:1:m_GF
                Estimated_BER(i,j) = min(max(LLR_PDF_Vec(i,j,1),1/(Num_Of_Iterations*(N/m))),1-1/(2^(m/m_GF)));
            end
        end
        
%         Estimated_Mutual_Information_Vec1 = sum(bsxfun(@times,((qfuncinv(max(Estimated_BER,1/(Num_Of_Iterations*(N/m))))).^2).',Estimated_Rates.'),1);
%         [~,Sorted_Estimated_SNR_Vec_Index_Vec1] = sort(Estimated_Mutual_Information_Vec1,'descend');
%         Sorted_By_Estimated_Constellation_Polar_Code_LLR1 = Constellation_Mapping_Arrays(Sorted_Estimated_SNR_Vec_Index_Vec1,:);
        
%         Estimated_Mutual_Information_Vec = zeros(1,size(Constellation_Mapping_Arrays,1));
%         Estimated_Mutual_Information_Vec1 = zeros(1,size(Constellation_Mapping_Arrays,1));
%         Estimated_Mutual_Information_Vec2 = zeros(1,size(Constellation_Mapping_Arrays,1));

%         for i=1:1:size(Constellation_Mapping_Arrays,1)
%             [d_min_square(i,:),Es(i,:)] = Calculate_Minimum_Distance(Constellation_Mapping_Arrays(i,:),m_BIPCM,m/m_GF);
%             d_min_square_to_Es(i,:) = d_min_square(i,:)./Es(i,:);
%             m_BIPCM_temp = m_BIPCM;
%             m_temp = m;
%             for j=1:1:m_GF
%                 if(m_BIPCM_temp>0)
%                     m_BIPCM_temp = m_BIPCM_temp - 1;
%                 else
%                     m_temp = m_temp - 1;
%                 end
%                 Estimated_Mutual_Information_Vec1(1,i) = Estimated_Mutual_Information_Vec(1,i)+()*log2(1+(qfuncinv(min(max(Estimated_BER(i,j),1/(Num_Of_Iterations*(N/m))),0.5)/(2^(m/m_GF)-1)).^2)/2);
%                 Estimated_Mutual_Information_Vec2(1,i) = Estimated_Mutual_Information_Vec(1,i)+0.5*log2(1+(qfuncinv(min(max(Estimated_BER(i,j),1/(Num_Of_Iterations*(N/m))),0.5)).^2)/2);
%                 Estimated_Mutual_Information_Vec(1,i) = Estimated_Mutual_Information_Vec(1,i)+0.5*log2(1+(qfuncinv(min(max(Estimated_BER(i,j),1/(Num_Of_Iterations*(N/m))),0.5)).^2)/(2*d_min_square_to_Es(i,j)*(log(2^m_temp)/log(2^(m/m_GF)))));
                
%                 Estimated_Mutual_Information_Vec(1,i) = Estimated_Mutual_Information_Vec(1,i) + log2(2^(m/m_GF)) - ((Estimated_BER(i,j)*log2(1/Estimated_BER(i,j))+(1-Estimated_BER(i,j))*log2(1/(1-Estimated_BER(i,j)))) + Estimated_BER(i,j)*log2(2^(m/m_GF)-1));
            
%             end
            
%             if(mod(i,Num_Of_Iterations_Print_State)==0)
%                 try
%                     if(and(is_sorted_by_N_m,is_sorted_by_Union_Bound))
%                         save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','Estimated_SER_Vec','Estimated_Is_Frozen_Bit_Index_Vec','LLR_Bin_Border_Vec','LLR_PDF_Vec','Estimated_BER');
%                     elseif(and(is_sorted_by_N_m,~is_sorted_by_Union_Bound))
%                         save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','LLR_Bin_Border_Vec','LLR_PDF_Vec','Estimated_BER');
%                     elseif(and(~is_sorted_by_N_m,is_sorted_by_Union_Bound))
%                         save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','Estimated_SER_Vec','Estimated_Is_Frozen_Bit_Index_Vec');
%                     end
%                 end
%                 sprintf('finish calculating label %d (out of %d labels)',i,size(Constellation_Mapping_Arrays,1))
%             end
            
%         end

        Estimated_Mutual_Information_Vec_Lower_Bound = (sum(log2(2^(m/m_GF)) - ((Estimated_BER.*log2(1./Estimated_BER)+(1-Estimated_BER).*log2(1./(1-Estimated_BER))) + Estimated_BER.*log2(2^(m/m_GF)-1)),2)).';
        
%         Estimated_Mutual_Information_Vec_Lower_Bound = (sum(log2(2^(m/m_GF)) - ((Estimated_BER.*log2(1./Estimated_BER)+(1-Estimated_BER).*log2(1./(1-Estimated_BER))) + Estimated_BER.*log2(2^(m/m_GF))),2)).';

        [~,Sorted_Estimated_SNR_Vec_Index_Vec_Lower_Bound] = sort(Estimated_Mutual_Information_Vec_Lower_Bound,'descend');
        Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound = Constellation_Mapping_Arrays(Sorted_Estimated_SNR_Vec_Index_Vec_Lower_Bound,:);
        
        Estimated_Mutual_Information_Vec_Upper_Bound = sum(0.5*log2(1+(qfuncinv(Estimated_BER)).^2/2),2).';

        [~,Sorted_Estimated_SNR_Vec_Index_Vec_Upper_Bound] = sort(Estimated_Mutual_Information_Vec_Upper_Bound,'descend');
        Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound = Constellation_Mapping_Arrays(Sorted_Estimated_SNR_Vec_Index_Vec_Upper_Bound,:);
        
    end

    try
        if(and(is_sorted_by_N_m,is_sorted_by_Union_Bound))
            save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','Estimated_BER','Union_Bound_Vec','Sorted_By_Union_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Estimated_Mutual_Information_Vec');
        elseif(is_sorted_by_N_m)
            save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','Estimated_BER','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Upper_Bound','Sorted_By_Estimated_Constellation_Polar_Code_LLR_Lower_Bound','Estimated_Mutual_Information_Vec');
        elseif(is_sorted_by_Union_Bound)
            save(fullfile(file_path,file_name),'Constellation_Mapping_Arrays','Estimated_BER','Union_Bound_Vec','Sorted_By_Union_Bound');
        end
    end
    
end